function processor = MSAprocessor
    processor.msa_processor=@msa_processor;
end

%% ========================== MSA PREPROCESSOR ===============================

function results=msa_processor(model)

    results = model;
    util = MSAutils;

    F = model.F;
    U = model.U;
    K = model.K;
    F_ext = model.F_ext;
    
    %Un = U;
    %Un(isnan(Un)) = 0;
    %F = K*Un;
    
    released = model.released;
    released(released==0)=[];


    idx = 1:numel(U);
    idx(released) = [];
    
    Un = U;
    Fn = F;
    Kn = K;
    T = eye(size(Kn));
%     if sum(released)>0
%         Kn(:,released) =[];
%         Kn(released,:) =[];
%         Fn(released) = [];
%         Un(released) = [];
%     end

    

    if ~strcmp(model.analysis_type,'beam')
        for i=1:model.number_supports
    
            support = model.support(i);
            t = support.angle;
            dof = support.node.dof;
        
            ids = support.node.directions;
    
            r = [cosd(t) -sind(t) 0;sind(t) cosd(t) 0;0 0 1];
            r = r(ids,ids);
    
            T(dof,dof)= r;
    
    
        end
    end

    T0 = T;
    Kr = T'*Kn*T;
    results.Kr = Kr;
    if sum(released)>0
        Kn(:,released) =[];
        Kn(released,:) =[];
        T(:,released) =[];
        T(released,:) =[];
        Fn(released) = [];
        Un(released) = [];
    end

    Krot = T'*Kn*T;
   
    %[Us,Fs] = solve_linear_system(Fn,Kn,Un);
    % TODO: Fn and Un is wrong when support is rotated...review it!!!
    [Urot,Frot] = solve_linear_system(Fn,Krot,Un);
    Us = T*Urot;
    Fs = T*Frot;
    
    F(idx) = Fs;
    U(idx) = Us;

    R = F - F_ext;

    results.T = T;
    
    Fr = F;
    Ur = U;
    Fr(idx) = Frot;
    Ur(idx) = Urot;
    
    Rr = Fr - T0'*F_ext;

    results.Ur = Ur;
    results.Fr = Fr;
    results.Rr = Rr;
    results.U = U;
    results.F = F;
    results.R = R;
    
    dofs = model.dof_node;
    ndof = numel(dofs);
    for i=1:model.number_elements
        
        %B = model.element(i).B;
        dof = model.element(i).dof;
        Ue = U(dof);
        Fe = F(dof);
        
        L(i) = results.element(i).L;
        if strcmp(model.analysis_type,'frame')
            Ke = results.element(i).K_global_fix;
            Fe(~isnan(Ue)) = nan(size(Fe(~isnan(Ue))));
       
            hinge = model.element(i).hinge;
    
            %Fh = Fe([2,3]);
            Fh = Fe([3,6]);
            Fh(hinge) = 0;
            Fe([3,6]) = Fh;
           
            Ue(~isnan(Fe)) = nan(size(Ue(~isnan(Fe))));
        
            [U_,F_] = solve_linear_system(Fe,Ke,Ue);
            Ulin(i,:) =  U_([1 2 4 5])';
            Uang(i,:) =  U_([3 6])';


        elseif strcmp(model.analysis_type,'truss')

            U_ = Ue;
            F_ = Fe;
            Ulin(i,:) =  U_([1 2 3 4])';
            Uang(i,:) =  Ulin(i,:);

        elseif strcmp(model.analysis_type,'beam')
            Ke = results.element(i).K_global_fix;
            Fe(~isnan(Ue)) = nan(size(Fe(~isnan(Ue))));
       
            hinge = model.element(i).hinge;
    
            Fh = Fe([2,4]);
            Fh(hinge) = 0;
            Fe([2,4]) = Fh;
           
            Ue(~isnan(Fe)) = nan(size(Ue(~isnan(Fe))));
        
            [U_,F_] = solve_linear_system(Fe,Ke,Ue);
            Ulin(i,:) =  U_([1 3])';
            Uang(i,:) =  U_([2 4])';

        elseif strcmp(model.analysis_type,'bar')
            U_ = Ue;
            F_ = Fe;
            Ulin(i,:) =  U_(1);
            Uang(i,:) =  Ulin(i,:);
        else

            U_ = Ue;
            F_ = Fe;

        end
          
        results.element(i).F = F_;
        results.element(i).U = U_;
        results.element(i).D = results.element(i).B*results.element(i).U;
        results.element(i).P = results.element(i).K_local*results.element(i).D;

       
        direc = {'x';'y';'xy'};
        P = results.element(i).P;
        Felem = results.element(i).F;
        D = results.element(i).D;
        Uelem = results.element(i).U;
        P1 = P(1:ndof);
        P2 = P(ndof+1:2*ndof);
        F1 = Felem(1:ndof);
        F2 = Felem(ndof+1:2*ndof);
        U1 = Uelem(1:ndof);
        U2 = Uelem(ndof+1:2*ndof);
        D1 = D(1:ndof);
        D2 = D(ndof+1:2*ndof);

        c = 1;
        for j=1:ndof
            %if id1(j)
   
            jj = dofs(j);
     
            %iforce(i,c,:) = msa2sm([P1(j) P2(j)],direc(j));
            results.element(i).force_local(j,:) = [P1(j) P2(j)];
            results.element(i).force_global(j,:) = [F1(j) F2(j)];
            results.element(i).displacment_local(j,:) = [D1(j) D2(j)];
            results.element(i).displacement_global(j,:) = [U1(j) U2(j)];
            results.element(i).internal_force(j,:) = util.msa2sm([P1(j) P2(j)],direc(jj));

            c = c + 1;
            %end
        end

    end

    if model.output_config.view.on_model{1}
         
        [opt, opt_type, coordinate_system] = util.get_direction(model);

        if ismember(opt_type,dofs)

            if strcmp(model.analysis_type,'beam')
                opt_type = opt_type - 1;
            end

            Qmin = zeros(model.number_elements,1);
            Qmax = zeros(model.number_elements,1);
            for i=1:model.number_elements
    
                if strcmp(coordinate_system,'local')
                    if strcmp(opt,'force')
    
                        Q = results.element(i).force_local(opt_type,:);
        
                    elseif strcmp(opt,'displacement')
    
                        Q = results.element(i).displacement_local(opt_type,:);
            
                    end
                elseif strcmp(coordinate_system,'global')
                    
                    if strcmp(opt,'force')
    
                        Q = results.element(i).force_global(opt_type,:);
        
                    elseif strcmp(opt,'displacement')

                        Q = results.element(i).displacement_global(opt_type,:);
            
                    end
    
                end
               
                results.element(i).Q = Q;
                Qmin(i) = min(Q);
                Qmax(i) = max(Q);
            end
       
        else
            error('MSAtool: cannot show the results in the select direction for this type of analysis.')
        end
        Qm = min(Qmin);
        QM = max(Qmax);
        %Flim_colorbar = linspace(Qm,QM,100);
        results.colorbar = [Qm QM];
    end

    if model.output_config.view.deformed

        l_max = max(L);
        u_max = max(abs(Ulin(:)));
        uang_max = max(abs(Uang(:)));
        
        % Update scale
        %scale = 0.1*l_max/u_max;
        scale = min(0.1*l_max/u_max,0.2/uang_max);
  
        for i=1:results.number_elements
            
            Be = results.element(i).B;
            B = eye(6);

            direc = [model.element(i).node_1.directions 3+model.element(i).node_2.directions]; 
     
            B(direc,direc)=Be;

            Uge = results.element(i).U;
            Ug = zeros(6,1);
            Ug(direc) = Uge;
      
            Ul = B*Ug;
            
            xc = results.element(i).node_1.x;
            yc = results.element(i).node_1.y;
            xf = results.element(i).node_2.x;
            yf = results.element(i).node_2.y;
            
            ds1 = Ul(3)*scale;
            ds2 = Ul(6)*scale;

            
            x1g = xc + Ug(1)*scale;
            x2g = xf + Ug(4)*scale;
            y1g = yc + Ug(2)*scale;
            y2g = yf + Ug(5)*scale;
    
            xgl = [x1g;y1g;ds1;x2g;y2g;ds2];
            gl = [x1g(1);y1g(1);0;x1g(1);y1g(1);0];
        
            xgl0 = xgl - gl;
            xloc0 = B*xgl0;
            
            x1 = xloc0(1);
            x2 = xloc0(4);
            y1 = xloc0(2);
            y2 = xloc0(5);
    
            xl = linspace(x1, x2, 100);

            if strcmp(model.analysis_type,'frame') || strcmp(model.analysis_type,'beam')
                A = [1 x1 x1^2 x1^3;0 1 2*x1 3*x1^2;1 x2 x2^2 x2^3;0 1 2*x2 3*x2^2];
                Y = [y1 ds1 y2 ds2]';
                C = A\Y;
                yl = (C(1) + C(2).*xl + C(3).*xl.^2 + C(4).*xl.^3);
                dyl = C(2) + 2*C(3).*xl + 3*C(4).*xl.^2;
            else
                A = [1 x1;1 x2];
                Y = [y1 y2]';
                C = A\Y;
                yl = (C(1) + C(2).*xl);
                dyl = ones(size(xl))*C(2);
            end
         
% 
%             B = variables(i).B;
%             gl = variables(i).gl;
%             xl = variables(i).xl;
%             yl = variables(i).yl;
%             dyl = variables(i).dyl;

            Bc = B(1:3,1:3);
            Z = Bc'*[xl;yl;dyl];
          
            Xg = Z + gl(1:3);
            xg = Xg(1,:);
            yg = Xg(2,:);
            
            results.element(i).x_deformed = xg;
            results.element(i).y_deformed = yg;
    
            results.node(results.element(i).node_1.id).x_deformed = xg(1);
            results.node(results.element(i).node_2.id).x_deformed = xg(end);
            results.node(results.element(i).node_1.id).y_deformed = yg(1);
            results.node(results.element(i).node_2.id).y_deformed = yg(end);
            
            % If view the data as a colormap on the deformed structure.
            if model.output_config.view.on_model{1}

                results.element(i).color_results = zeros(100,3);

                Qe = results.element(i).Q;
                Qmin_e = Qe(1);
                Qmax_e = Qe(2);

                Fv = linspace(Qmin_e,Qmax_e,100);
                if strcmp(coordinate_system,'global')
                    Flim = linspace(Qm,QM,100);
                else
                    Flim = Fv;
                end

                results.element(i).color_results_limits = Flim;
                
                if all(Flim(2:end) ~= Flim(1))
                    
      
                    c0 = jet(100);
                    for ii=1:100
                        results.element(i).color_results(ii,:) = util.interpcolor(Flim,c0,Fv(ii));
                    end

                else
                    
                    for ii=1:100
                        results.element(i).color_results(ii,:) = [0 1 0];
                    end
                end

            end


        end

    end
    
end

function [U,F] = solve_linear_system(F,K,U)

    idx1 = isnan(F);
    idx2 = isnan(U);
    
    rows = 1:size(K,1);
    cols = 1:size(K,2);
    Fp = F(idx1==0);
    us = U(idx2==0);
    
    rows_p = rows(idx1==0);
    cols_s = cols(idx2==0);
    rows_s = rows(idx1==1);
    cols_p = cols(idx2==1);
    
    
    Kpp = K(rows_p,cols_p);
    Kps = K(rows_p,cols_s);
    Ksp = K(rows_s,cols_p);
    Kss = K(rows_s,cols_s);
    up = Kpp\(Fp - Kps*us);
    Fs = Ksp*up + Kss*us;
    
    U(idx2)=up;
    F(idx1)=Fs;


end

