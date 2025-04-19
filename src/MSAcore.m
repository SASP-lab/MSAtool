function core = MSAcore
    core.msa_read_data=@msa_read_data;
    core.msa_processor=@msa_processor;
    core.msa_results=@msa_results;
end

%% ========================== MSA READ DATA ===============================

function model=msa_read_data(varargin)

    % University of Minnesota Twin Cities
    % Department of Civil, Environmental, and Geo- Engineering
    % Author: Prof. Ketson R. M. dos Santos
    % Email: dossantk@umn.edu
    % Date of last modification: Mar 13, 2024.
    % Version 1.0.0
    %
    % READ_DATA.m: 
    % Read the data from the input file to perform the analysis of framed
    % structures.

    util = MSAutils;

    %set functions
    test_line_size = util.test_line_size;
    model = {};
    
    if numel(varargin)==1
        input_file=varargin{1};
        verbose = true;
    elseif numel(varargin) == 2
        input_file=varargin{1};
        verbose=varargin{2};
    else
        error('Error: the number of input variables in read_data is larger than expected')
    end
    
    model.input_file = input_file;
    if exist(input_file)
         if verbose
              util.header(1,model.input_file,[],[]);
         end
    else
         % File does not exist.
         error('Error: input file does not exist!');
    end

    file = input_file;
    fid = fopen(file,'r');
    
    %lines_input = textread(file,'%s','delimiter','\r');
    lines_scan = textscan(fid,'%s','delimiter','\r');
    lines_input = lines_scan{1};

    F = [];
    U = [];
    Pd = [];
    node = {};
    element = {};
    support = {};
    num_nodes = nan;
    num_elements = nan;

    output_config = {};
    output_config.save_output = false;
    output_config.show_output = false;

    numlines = length(lines_input);
    iline = 1;
    while iline<numlines
    
        line = lines_input{iline};
        if strcmp(line,'*ANALYSIS')
            analysis_type = lines_input{iline+1};

            if verbose
                disp(['* Analysis type: [' 8 [analysis_type ']'] 8])
            end

            iline=iline+1;
        end

        if strcmp(analysis_type,'bar_1D')

            num_dof_node = 1;
            
            xyz = [1 0 0];
            node_dof = [1 0 0 0 0 0];
            [iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
                       read_element(iline,line,lines_input,node,element,...
                     support,num_nodes,num_elements,F,Pd,U,node_dof,analysis_type,util,xyz);

            %error('Error: bar not implemented!')

        elseif strcmp(analysis_type,'beam_2D')

            num_dof_node = 2;

            xyz = [1 0 0];
            node_dof = [0 1 0 0 0 1];
            [iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
                       read_element(iline,line,lines_input,node,element,...
                     support,num_nodes,num_elements,F,Pd,U,node_dof,analysis_type,util,xyz);


            %error('Error: beam not implemented!')

        elseif strcmp(analysis_type,'truss_2D')

            num_dof_node = 2;

            xyz = [1 1 0];
            node_dof = [1 1 0 0 0 0];
            [iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
                       read_element(iline,line,lines_input,node,element,...
                     support,num_nodes,num_elements,F,Pd,U,node_dof,analysis_type,util,xyz);

            
            %error('Error: truss_2D not implemented!')

        elseif strcmp(analysis_type,'frame_2D')
    
             num_dof_node = 3;

             xyz = [1 1 0];
             node_dof = [1 1 0 0 0 1];
            [iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
                       read_element(iline,line,lines_input,node,element,...
                     support,num_nodes,num_elements,F,Pd,U,node_dof,analysis_type,util,xyz);

            %error('Error: frame_2D not implemented!')

         elseif strcmp(analysis_type,'grid_3D')
    
             num_dof_node = 3;

             xyz = [1 0 1];
             node_dof = [0 1 0 1 0 1];
            %[iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
            %           read_element(iline,line,lines_input,node,element,...
            %         support,num_nodes,num_elements,F,Pd,U,node_dof,analysis_type,util,xyz);

            error('Error: frame_2D not implemented!')

        elseif strcmp(analysis_type,'frame_3D')

            error('Error: frame_3D not implemented!')


        else

            error('Error: analysis not recognized!')

        end

        show_output = false; % default.
        
        if strcmp(line,'*SHOW_OUTPUT')
            show_output = str2double(lines_input{iline+1});
            
            if show_output ~= 0 && show_output ~= 1
                error('Error: show output must be either 0 or 1.')
            else
                output_config.show_output = logical(show_output);
            end
            
            iline=iline+1;
        end
        
        save_output = false; % default.
        %output_config.save_output = save_output;
        if strcmp(line,'*SAVE_OUTPUT')

            save_output = str2double(lines_input{iline+1});
            
            if save_output ~= 0 && save_output ~= 1
                error('Error: save output must be either 0 or 1.')
            else
                output_config.save_output = logical(save_output);
            end
            
            iline=iline+1;
        end
    
        if strcmp(line,'*VISUALIZATION')
            visual_bool = logical(str2double(lines_input{iline+1}));
            
            if verbose
                str_fig = 'FALSE';
                if visual_bool
                    str_fig = 'TRUE';
                end

                disp(['* Show model and results: [' 8 [str_fig ']'] 8])
            end

            if numel(visual_bool) ~=1
                error('Error in *VISUALIZATION')
            end

            output_config.visual_bool = visual_bool;
            
            if visual_bool

                visual_in = lines_input{iline+2};
                visual_ = split(visual_in);
           
                bool_fig = logical(str2double(visual_{1}));
                output_config.fig_results_bool = bool_fig;
                
                visual = visual_(2:end);
               
                test_line_size(visual_, 4, '*VISUALIZATION');
    
                info = split(lines_input{iline+3});
                test_line_size(info, 4, '*VISUALIZATION');
                output_config.show_id = logical(str2double(info{1}));
                output_config.show_forces = logical(str2double(info{2}));
                output_config.show_reactions = logical(str2double(info{3}));
                output_config.show_deformed = logical(str2double(info{4}));
                
                elem_info = split(lines_input{iline+4});
                
                element_id = zeros(numel(elem_info), 1);
                for i=1:numel(elem_info)
                    element_id(i) = str2double(elem_info(i));
                end
                element_id(isnan(element_id))=[];
                output_config.element_id = element_id;

                if bool_fig

                    if strcmp(visual{1},'displacements') || strcmp(visual{1},'forces')
                        
                        output_config.results_on_fig = visual{1};
        
                        if strcmp(visual{2},'x')
        
                            output_config.force_dir = 1;
        
                        elseif strcmp(visual{2},'y')
        
                            output_config.force_dir = 2;
        
                        elseif strcmp(visual{2},'z')
        
                            output_config.force_dir = 3;
        
                        elseif strcmp(visual{2},'zy')
        
                            output_config.force_dir = 4;
        
                        elseif strcmp(visual{2},'zx')
        
                            output_config.force_dir = 5;
        
                        elseif strcmp(visual{2},'xy')
        
                            output_config.force_dir = 6;
                        else
                            error('Error: inconsistent information in *VISUALIZATION (direction).');
        
                        end
                       
                        if strcmp(visual{3},'local') || strcmp(visual{3},'global')
                            output_config.coordinate_system = visual{3};
                        else
                            error('Error: inconsistent information in *VISUALIZATION (global/local).');
                        end
        
                    else
                        output_config.force_dir = NaN;
                        output_config.results_on_fig = NaN;
                    end
                    
                    if verbose
                        disp(['* Results: [' 8 [visual_in ']'] 8])
                    end
                else
                    output_config.force_dir = NaN;
                    output_config.results_on_fig = NaN;
                    %output_config.show_id = NaN;
                    %output_config.show_forces = NaN;
                    output_config.coordinate_system = NaN;
                     output_config.results_on_fig = false;
                end
    
                iline=iline+4;

            else
                iline=iline+1;
            end
            
        end
      
        if strcmp(line,'*END')
            disp(line);
            break
        end
        iline=iline+1;
    end

    num_elements = numel(element);
    num_nodes = numel(node);
    model.analysis_type = analysis_type;
    model.number_elements = num_elements;
    model.number_nodes = num_nodes;
    model.element = element;
    model.node = node;
    model.support = support;
    model.number_supports = numel(support);
    model.output_config = output_config;

    errors = test_model(model);
    if ~isempty(errors)

        fprintf('Input data is inconsistent!');
        for i=1:numel(errors)
            fprintf('%s\n',errors(i));
        end
    end
    

    if visual_bool
        if model.output_config.fig_results_bool
            if sum(sum(element_id==1:num_elements,2))~=numel(element_id)
                error('Error: element ID in *VISUALIZATION.');
            end
        end
    end
   
    K = get_global_system(element, num_elements, num_nodes, num_dof_node);
    model.K = K;
    model.F = F;
    model.U = U;
    model.Pd = Pd;

    if verbose
        disp('   ');
    end
    
end

function [iline,node,element,F,Pd,U,num_nodes,num_elements,support]=...
    read_element(iline,line,lines_input,node,element,...
                    support,num_nodes,num_elements,F,Pd,U,node_dof,analysis,util,xyz)

    test_line_size = util.test_line_size;
    
    ndof = sum(node_dof);
    dim = str2num(analysis(end-1:end-1));

    if strcmp(line,'*NODES')
        
        num_nodes = str2double(lines_input{iline+1});
        node = struct.empty(num_nodes,0);
        loc_id_vec  = zeros(num_nodes,1);
        iline=iline+1;
        
        F = zeros(ndof*num_nodes, 1);
        Pd = zeros(ndof*num_nodes, 1);
        % count = 1;
        for i=1:num_nodes

            node_data = split(lines_input{iline+1});
        
            %test_line_size(node_data, 1+dim+2*ndof, '*NODES');
            loc_id = str2double(node_data{1});
            
            if any(loc_id_vec(:) == loc_id) || loc_id<=0
                error('Error: node id inconsistent in *NODE.')
            else
                loc_id_vec(i) = loc_id;
            end

            node(loc_id).id = str2double(node_data{1});
            node(loc_id).node_dof = node_dof;
            
            position = zeros(3,1);

            cc=2;
            for j=1:dim
                position(j) = str2double(node_data{cc});
                cc = cc + 1;
            end
    
            force = zeros(6,1);
            displacement = zeros(6,1);
            for j=1:6
                if node_dof(j)==1
                    force(j) = str2double(node_data{cc});
                    displacement(j) = str2double(node_data{cc+ndof});
                    cc = cc + 1;
                end
            end

            node(loc_id).x = position(1);
            node(loc_id).y = position(2);
            node(loc_id).z = position(3);
            

            % Nodal prescribed displacements.
            node(loc_id).d = displacement;

            % Nodal forces
            node(loc_id).f = force;
            node(loc_id).hinge = 0;
            
            force_ = force(node_dof==1);
            displacement_ = displacement(node_dof==1);
            loc_f = (loc_id*ndof-ndof+1):(loc_id*ndof);

            F(loc_f) = force_;
            Pd(loc_f) = displacement_;
            
            node(loc_id).number_dof = ndof;
            node(loc_id).D_node = zeros(ndof,6);
            
            for j=1:ndof
                node(loc_id).dof_global(j) = ndof*node(loc_id).id - (ndof-j);
                node(loc_id).D_node(j,node_dof==1) = 1;
            end

            iline=iline+1;

        end

    end

    if strcmp(line,'*SUPPORTS')

        num_supports = str2double(lines_input{iline+1});
        if num_supports > num_nodes
            error('Error: number of supports is larger than the number of nodes!')
        end

        support = struct.empty(num_supports,0);
        iline=iline+1;

        U = nan(ndof*num_nodes,1);
        for i=1:num_supports
            support_data = split(lines_input{iline+1});
            test_line_size(support_data, 1+ndof, '*SUPPORTS');
            
            support(i).node_id = str2double(support_data{1});
            constraint = ones(1,6);
    
            cc = 2;
            for j=1:6
                if node_dof(j)==1
                    constraint(j) = str2double(support_data{cc}); 
                    cc = cc + 1;
                end
            end
            support(i).angle = 0;
          
            idx = 1:6;
            const_local = constraint(idx(node_dof==1));
            node_id = support(i).node_id;
            %dof_node = node(node_id).dof_global;
            support(i).node = node(node_id);
            dof_global = node(node_id).dof_global;
            %support(i).dof_global_constraint = dof_global(1);
            support(i).dof_global_constraint = dof_global(logical(const_local));
            
            %F(dof_node) = nan;
            %U(dof_node) = 0; 

            F(support(i).dof_global_constraint) = nan;
            U(support(i).dof_global_constraint) = 0; % no support displacement. Future implementation...

            support(i).constraint = constraint;
            %support(i).constraint(2) = 1;
            %support(i).constraint(3) = 1;

            iline=iline+1;
        end

    end

    if strcmp(line,'*ELEMENTS')
        num_elements = str2double(lines_input{iline+1});
    
        element = struct.empty(num_elements,0);
      
        iline=iline+1;
        for i=1:num_elements
            element_data = split(lines_input{iline+1});

            test_line_size(element_data, 6, '*ELEMENTS');
            
            element(i).id = str2double(element_data{1});
            element(i).E = str2double(element_data{2});
            element(i).A = str2double(element_data{3});
            element(i).I = str2double(element_data{4});
           
            id_1 = str2double(element_data{5});
            id_2 = str2double(element_data{6});
         
            element(i).node_1 = node(id_1);
            element(i).node_2 = node(id_2); 

            if dim==1
                dx = element(i).node_2.x - element(i).node_1.x;
                element(i).L = abs(dx);
                element(i).angle = 0;
            elseif dim==2
                dx = element(i).node_2.x - element(i).node_1.x;
                dy = element(i).node_2.y - element(i).node_1.y;
                element(i).angle = atan2d(dy,dx);
                element(i).L = sqrt(dx^2 + dy^2);
            elseif dim==3
                error('Error: 3D analysis not implemented yet, wait for the next version.')
            else
                error('Error: analysis must be 1D, 2D, or 3D.')
            end

            if strcmp(analysis,'bar_1D')
                [K_local,D,dof_vec] = get_matrices_bar(element(i), num_nodes);
                K_global = K_local;
                B = eye(2);
            elseif strcmp(analysis,'beam_2D')    
                [K_local,D,dof_vec] = get_matrices_beam(element(i), num_nodes);
                K_global = K_local;
                B = eye(4);
            elseif strcmp(analysis,'truss_2D')
                [K_local,K_global,B,D,dof_vec]  = get_matrices_truss_2d(element(i), num_nodes);
            elseif strcmp(analysis,'frame_2D')
                [K_local,K_global,B,D,dof_vec] = get_matrices_frame_2d(element(i), num_nodes);
            else
                error('Error: analysis not implemented yet, wait for the next version.')
            end

            element(i).K_local = K_local;
            element(i).K_global = K_global;
            element(i).B = B;
            element(i).D = D;
            element(i).dof_global = dof_vec;
            
            iline=iline+1;
        end
    end

end

function errors = test_model(model)

    analysis = model.analysis_type;
    errors = [];
    if strcmp(analysis,'bar_1D')
        
        if model.number_nodes ~= (model.number_elements + 1)
            errors = [errors 'Model verification: number of nodes and elements are inconsistent!'];
        end

    elseif strcmp(analysis,'beam_2D')

    elseif strcmp(analysis,'truss_2D')

    elseif strcmp(analysis,'frame_2D')

    else
        errors = [errors 'Model verification: *ANALYSIS not recognized/implemented!'];
    end

end

function [K_local,D,dof_vec] = get_matrices_bar(element, num_nodes)

    EA = element.E*element.A;
    L = element.L;
    K_local = [EA/L -EA/L;
               -EA/L EA/L];
    
    dof_vec = [element.node_1.dof_global element.node_2.dof_global];
    D = zeros(2, num_nodes);

    for i=1:2
        D(i,dof_vec(i)) = 1;
    end

end

function [K_local,D,dof_vec] = get_matrices_beam(element, num_nodes)

    EI = element.E*element.I;
    L = element.L;
    K_local = [12*EI/L^3 6*EI/L^2 -12*EI/L^3 6*EI/L^2;
               6*EI/L^2 4*EI/L -6*EI/L^2 2*EI/L;
               -12*EI/L^3 -6*EI/L^2 12*EI/L^3 -6*EI/L^2;
               6*EI/L^2 2*EI/L -6*EI/L^2 4*EI/L];

    dof_vec = [element.node_1.dof_global element.node_2.dof_global];
    D = zeros(4, 2*num_nodes);
    
    for i=1:4
        D(i,dof_vec(i))=1;
    end
    

end

function [K_local,K_global,B,D,dof_vec] = get_matrices_truss_2d(element, num_nodes)
    

    EA = element.E*element.A;
    L = element.L;
    K_local = [EA/L 0 -EA/L 0;
               0 0 0 0;
               -EA/L 0 EA/L 0;
               0 0 0 0];

    c = cosd(element.angle);
    s = sind(element.angle);
    B = [c s 0 0;
         -s c 0 0;
         0 0 c s;
         0 0 -s c];

    K_global = B'*K_local*B;

    dof_vec = [element.node_1.dof_global element.node_2.dof_global];
    D = zeros(4, 2*num_nodes);
    
    for i=1:4
        D(i,dof_vec(i))=1;
    end
   
end

function [K_local,K_global,B,D,dof_vec] = get_matrices_frame_2d(element, num_nodes)

    EI = element.E*element.I;
    EA = element.E*element.A;
    L = element.L;
    K_local = [EA/L 0 0 -EA/L 0 0;
        0 12*EI/L^3 6*EI/L^2 0 -12*EI/L^3 6*EI/L^2;
        0 6*EI/L^2 4*EI/L 0 -6*EI/L^2 2*EI/L;
        -EA/L 0 0 EA/L 0 0;
        0 -12*EI/L^3 -6*EI/L^2 0 12*EI/L^3 -6*EI/L^2;
        0 6*EI/L^2 2*EI/L 0 -6*EI/L^2 4*EI/L];

    c = cosd(element.angle);
    s = sind(element.angle);
    B = [c s 0 0 0 0;
         -s c 0 0 0 0;
          0 0 1 0 0 0;
          0 0 0 c s 0;
          0 0 0 -s c 0;
          0 0 0 0 0 1];

    K_global = B'*K_local*B;
    
    dof_vec = [element.node_1.dof_global element.node_2.dof_global];
    D = zeros(6, 3*num_nodes);
    
    for i=1:6
        D(i,dof_vec(i))=1;
    end

end

function [K_local,K_global,B,D,dof_vec] = get_matrices_grid_2d(element, num_nodes)

    EI = element.E*element.I;
    GJ = element.G*element.J;
    L = element.L;
    K_local = [GJ/L 0 0 -GJ/L 0 0;
        0 12*EI/L^3 6*EI/L^2 0 -12*EI/L^3 6*EI/L^2;
        0 6*EI/L^2 4*EI/L 0 -6*EI/L^2 2*EI/L;
        -GJ/L 0 0 GJ/L 0 0;
        0 -12*EI/L^3 -6*EI/L^2 0 12*EI/L^3 -6*EI/L^2;
        0 6*EI/L^2 2*EI/L 0 -6*EI/L^2 4*EI/L];

    c = cosd(element.angle);
    s = sind(element.angle);
    B = [c 0 s 0 0 0;
         0 1 0 0 0 0;
        -s 0 c 0 0 0;
          0 0 0 c 0 s;
          0 0 0 0 1 0;
          0 0 0 -s 0 c];

    K_global = B'*K_local*B;
    
    dof_vec = [element.node_1.dof_global element.node_2.dof_global];
    D = zeros(6, 3*num_nodes);
    
    for i=1:6
        D(i,dof_vec(i))=1;
    end

end

function K = get_global_system(element,num_elements, num_nodes, num_dof_node)

    K = zeros(num_dof_node*num_nodes, num_dof_node*num_nodes);
    for i=1:num_elements
        K = K + element(i).D'*element(i).K_global*element(i).D;
    end

end

%% ========================= MSA PROCESSOR ================================

function results = msa_processor(model)
   
    % Read the utility functions.
    util = MSAutils;
    sort_back = util.sort_back;

    % Read the option to save output.
    save_output = model.output_config.save_output;

    % Create the structure for store the results data.
    results = {};
    results.save_output = save_output;

    % Determine the analysis type.
    analysis_type = model.analysis_type;
    results.input_file = model.input_file;

    if save_output

        getdate = string(datetime('now','Format','MMMM-dd-yyyy_HH_mm_ss'));
        output_file = strcat(analysis_type,'_output_',getdate,'.txt');
        fileID = fopen(output_file,'w');
        util.header(fileID,model.input_file,analysis_type,output_file);
        results.output_file = output_file;
        
    else
        fileID = 0;
    end

    F = model.F;
    U = model.U;
    K = model.K;
    Pd = model.Pd;
    
    % Check this implementation.
    U(Pd~=0)=Pd(Pd~=0);
    F0 = F(Pd~=0);
    F(Pd~=0)=NaN;
    Fpd = K*Pd;
    %Fpd(Pd==0)=0;
   
    %F = F - Fpd;
   
    results.K = K;
    
    % Sort the vector of force to separate known and unknown quantites.
    [F_sort, F_idx] = sort(F,'ascend');
    U_sort = U(F_idx);
    K_sort = K(F_idx,F_idx);
    
    idx = 1:numel(U);
    id_p = idx(~isnan(F_sort));
    id_s = idx(isnan(F_sort));
    
    Kpp = K_sort(id_p, id_p);
    Kps = K_sort(id_p, id_s);
    Ksp = K_sort(id_s, id_p);
    Kss = K_sort(id_s, id_s);
    us = U_sort(id_s);
    Fp = F_sort(id_p);
    
    % Find displacements:
    up = Kpp\(Fp - Kps*us);
    Fs = Ksp*up + Kss*us;
  
    F = cat(1,Fp,Fs); 
    F = sort_back(F, F_idx, 1 );
    U = cat(1,up,us);
    U = sort_back(U, F_idx, 1 );

    %U = U + Pd;
    %F = F + Fpd;

    results.element = model.element;
    results.node = model.node;
    results.support = model.support;
    results.number_nodes = model.number_nodes;
    results.number_supports = model.number_supports;
    results.number_elements = model.number_elements;
    results.analysis_type = model.analysis_type;
   
    for i=1:model.number_elements

        %ndof = model.element(i).node_1.number_dof;
        results.element(i).U = U(model.element(i).dof_global);
        delta = model.element(i).B*U(model.element(i).dof_global);
        results.element(i).delta = delta;
        results.element(i).P = model.element(i).K_local*delta;
        results.element(i).F = model.element(i).K_global*U(model.element(i).dof_global);
        %results.element(i).P
        %model.element(i).dof_global

    end
    
    ndof = results.node(1).number_dof;
    soma = zeros(ndof,1);

    if sum(abs(soma))==0
        fprintf('* Equilibrium: OK!\n');
    else
        error('Error: equilibrium not satisfied!')
    end

    c = 1;
    for i=1:model.number_nodes
        for j=0:ndof-1
            soma(j+1) = soma(j+1) + F(c+j);
        end
        c = c + ndof;
    end

    %F = F + Fpd;
    results.F = F;
    results.U = U;
    
    %v = 1:6;
    for i=1:model.number_supports

        id_constr = model.support(i).constraint.*model.support(i).node.node_dof;

        Q = model.support(i).node.f;
        Q(isnan(Q)) = 0;
        Q(id_constr == 0) = 0;
     
        id_const = model.support(i).dof_global_constraint;
        node_ndof = model.support(i).node.number_dof;
        node_dof_global = model.support(i).node.dof_global;
        R0 = zeros(1,node_ndof);
        
        [~,posi]=ismember(id_const,node_dof_global);
        posi(posi==0)=[];
        
        %R0(id_const == node_dof_global) = F(id_const);
        R0(posi) = F(id_const);
       
        R1 = [0 0 0 0 0 0];
        R1(model.support(i).node.node_dof==1) = R0;

        results.support(i).R = R1 - Q';
       
        %results.support(i).R = F(id_const);
        %results.support(i).U = U(id_const);
    end

    %p_reac = [];
    %for i=1:model.number_supports
    %    dof=model.support(i).node.dof;
        %id_dof=model.support(i).constraint([1,2]); %Truss only
        %p_reac = [p_reac dof(logical(id_dof))];
    %end
    
    % show_save
    util.show_save_results(results,model,fileID)
  
end

%% ========================== MSA RESULTS =================================

function msa_results(model,results)
    
    % msa_resuls: show the results obtained after solving the linear
    % analysis.

    % Get informationa from model.
    show_id = model.output_config.show_id;
    show_forces = model.output_config.show_forces;
    fig_results_bool = model.output_config.fig_results_bool;
    
    show_reactions = model.output_config.show_reactions;
    show_deformed = model.output_config.show_deformed;
    results_on_fig = model.output_config.results_on_fig;
    force_dir = model.output_config.force_dir;
    elem_id = model.output_config.element_id;
    coordinate_system = model.output_config.coordinate_system;

    % Call the utility functions.
    util = MSAutils;
    %get_force_internal = util.get_force_internal;
    interpcolor = util.interpcolor;

    % logical variable to show results as a colormap on figures.
    show_on_fig = strcmp(results_on_fig,'forces') || strcmp(results_on_fig,'displacements');
    
    % Determine the type of analysis and the dimension of the problem.
    analysis = model.analysis_type;
    dim = str2double(analysis(end-1:end-1));
    
    % Get the nodes' coordinates.
    x = zeros(1,model.number_nodes);
    y = zeros(1,model.number_nodes);
    z = zeros(1,model.number_nodes);
    for i=1:model.number_nodes
        x(i) = model.node(i).x;
        y(i) = model.node(i).y;
        z(i) = model.node(i).z;
    end
    
    % Find the minimum and maximum coordinates.
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    minz = min(z);
    maxz = max(z);
    
    mxx = max([maxx,maxy,maxz]);
    mix = max([minx,miny,minz]);

    maxx = mxx;
    minx = mix;
    maxy = mxx;
    miny = mix;
    maxz = mxx;
    minz = mix;

    % Determine the size of each dimension.
    Dx = maxx - minx;
    Dy = maxy - miny;
    Dz = maxz - minz;
    Dxyz = [Dx Dy Dz];

    % Create a box to plot the data.
    box = [minx maxx miny maxy minz maxz];

    % Create logical variables to determine if you must show the results on
    % figure.
    if show_id
        if show_on_fig
            node_num_1 = false;
            elem_num_1 = false;
            node_num = true;
            elem_num = true;
        else
            node_num_1 = true;
            elem_num_1 = true;
            node_num = false;
            elem_num = false;
        end

    else
        node_num_1 = false;
        elem_num_1 = false;
        node_num = false;
        elem_num = false;
    end




    % If show on figure, plot the model with high transparency, otherwise
    % show the elements with black colors.
    if show_on_fig
        show_model(model,show_forces,[0.8,0.8,0.8],'k','k',false,box,node_num_1,elem_num_1,0.1);
    else
        show_model(model,show_forces,[0,0,0],'k','k',false,box,node_num_1,elem_num_1,1);
    end
    
    % If show the deformed shape of the structure.
    if show_deformed

        % Get the nodal displacements.
        U = results.U;
        
        Uval = [];
        for i=1:results.number_nodes

            Ui = U(results.node(i).dof_global);
            if dim==1
                Uval = [Uval Ui(1)];
            elseif dim==2
                Uval = [Uval Ui(1:2)'];
            else
                %Uval = [Uval Ui(1:3)'];
            end

        end

        Umax = max(abs(Uval));
        if Umax == 0
            scale = 0;
        else
            scale = Dx/(15*Umax);
        end

        U = scale*U;
        for i=1:results.number_nodes
            
            Ui = U(results.node(i).dof_global);
            
            if dim==1
                
                results.node(i).x = results.node(i).x + Ui(1);

             elseif dim==2
                 if strcmp(analysis,'beam_2D')
                    results.node(i).y = results.node(i).y + Ui(1);
                 else
                    results.node(i).x = results.node(i).x + Ui(1);
                    results.node(i).y = results.node(i).y + Ui(2);
                 end
 
             else
 
                 %results.node(i).x = results.node(i).x + Ui(1)*scale;
                 %results.node(i).y = results.node(i).y + Ui(2)*scale;
                 %results.node(i).z = results.node(i).z + Ui(3)*scale;

            end

        end
      
        for i=1:results.number_elements
            results.element(i).node_1 = results.node(results.element(i).node_1.id);
            results.element(i).node_2 = results.node(results.element(i).node_2.id);
        end

        for i=1:numel(results.support)
            node_id = results.support(i).node_id;
            results.support(i).node = results.node(node_id);
        end
        
        if show_on_fig
            show_model(results,false,[0,0,0],'g','g',true,box,node_num,elem_num,0);
        end
        
    end

    if show_reactions

        for i=1:results.number_supports
     
            force = results.support(i).R;
            node = results.support(i).node;
            position = [node.x, node.y, node.z];
            plot_forces(force,position,Dxyz,dim,'r');

        end
    
    end

    if show_on_fig

        FM = zeros(results.number_elements,1);
        Fm = zeros(results.number_elements,1);
        
        Ft = [];
        for i=1:results.number_elements

            ndof = results.element(i).node_1.number_dof;
            node_dof = results.element(i).node_1.node_dof; 
            node_dof(node_dof==1)=1:ndof;
          
            force_dir_ = node_dof(force_dir);
            
            if strcmp(coordinate_system,'local')

                if strcmp(results_on_fig,'forces')
                    
                    %Pm = get_force_internal(results.element(i));
                    Pm = results.element(i).P;

                    Fo = Pm([force_dir_ ndof + force_dir_]);
    
                elseif strcmp(results_on_fig,'displacements')
    
                    Pm = results.element(i).delta;
                    Fo = Pm([force_dir_ ndof + force_dir_]);
    
                end

            elseif strcmp(coordinate_system,'global')

                if strcmp(results_on_fig,'forces')
                    
                    Pm = results.element(i).F;
                    Fo = Pm([force_dir_ ndof + force_dir_]);
    
                elseif strcmp(results_on_fig,'displacements')
    
                    Pm = results.element(i).U;
                    Fo = Pm([force_dir_ ndof + force_dir_]);
    
                end

            end
        
            FM(i) = max(Fo);
            Fm(i) = min(Fo);

            ndof = numel(results.element(i).dof_global);
            t1 = results.element(i).U(ndof/2);
            t2 = results.element(i).U(ndof);
            Ft = [t1 t2 Ft];

        end
    
        Tmax = max(abs(Ft));
        Fmax = max(FM);
        Fmin = min(Fm);
        Flim = linspace(Fmin,Fmax,100);
      
        if Tmax~=0
            scale_rot=0.3/tan(Tmax);
        else
            scale_rot=0.3;
        end
    
        %scale_rot = scale;
        %scale_rot = 1;
     
        elems = {};
        for i=1:results.number_elements
         
            ndof = results.element(i).node_1.number_dof;

            x10 = model.element(i).node_1.x;
            y10 = model.element(i).node_1.y;
            z10 = model.element(i).node_1.z;
            x20 = model.element(i).node_2.x;
            y20 = model.element(i).node_2.y;
            z20 = model.element(i).node_2.z;

            x1 = results.element(i).node_1.x;
            y1 = results.element(i).node_1.y;
            z1 = results.element(i).node_1.z;
            x2 = results.element(i).node_2.x;
            y2 = results.element(i).node_2.y;
            z2 = results.element(i).node_2.z;

            dx1 = x1 - x10; 
            dy1 = y1 - y10; 
            dz1 = z1 - z10; 
            dx2 = x2 - x20; 
            dy2 = y2 - y20; 
            dz2 = z2 - z20; 

            if strcmp(coordinate_system,'local')
                if strcmp(results_on_fig,'forces')
                    %Pm = get_force_internal(results.element(i));
                    Pm = results.element(i).P;
                elseif strcmp(results_on_fig,'displacements')
                    Pm = results.element(i).delta;
                end
            elseif strcmp(coordinate_system,'global')
                if strcmp(results_on_fig,'forces')
                    Pm = results.element(i).F;
                elseif strcmp(results_on_fig,'displacements')
                    Pm = results.element(i).U;
                end
            end
    
            elems(i).Pm = Pm;
            
            if force_dir > 6 || force_dir < 0
                error('Error: error in force_dir!');
            end
    
            Pdir = Pm([force_dir_ ndof + force_dir_]);
            elems(i).Pdir = Pdir;
            
            if all(Flim(2:end) ~= Flim(1))

                if dim<=2
                    
                    %scale=0.2/Tmax;
                    element = results.element(i);
                    [x,y]=util.coordinates_deformed(element,x10,y10,x1,x2,y1,y2,dx1,dx2,dy1,dy2,scale_rot);
                    
                    c0 = jet(100);
                    Fv = linspace(Pdir(1),Pdir(2),100);
                    c = zeros(100,3);
                    for ii=1:100
                        c(ii,:)=interpcolor(c0,Flim,Fv(ii));
                    end
                    
                else
                    disp('Results on figure not inmplemented for 3D structures!')
                end

            else

                if dim<=2

                    %scale=0.2/Tmax;
                    element = results.element(i);
                    [x,y]=util.coordinates_deformed(element,x10,y10,x1,x2,y1,y2,dx1,dx2,dy1,dy2,scale_rot);

                    c = zeros(100,3);
                    for ii=1:100
                        c(ii,:)=[0 1 0];
                    end
                
                else
                    disp('Results on figure not inmplemented for 3D structures!')
                end
            end
        
            n = length(x);
            fv = [1:n-1;2:n]';
            patch('faces',fv,'vertices',[x; y]',...
                  'faceVertexCData',c,...
                  'edgecolor','flat',...
                  'linewidth',2);
            
            colorbar;
            colormap jet;
            if Fmin ~= Fmax    
                set(gca,'CLim',[Fmin Fmax]);
            else
                if Fmin==0
                    set(gca,'CLim',[-1 1]);
                else
                    set(gca,'CLim',[-Fmin Fmin]);
                end
            end
    
        end


        if numel(elem_id)>0
    
            if strcmp(results_on_fig,'forces')
                plots(1).title = 'Force x';
                plots(2).title = 'Force y';
                plots(3).title = 'Force z';
                plots(4).title = 'Moment x';
                plots(5).title = 'Moment y';
                plots(6).title = 'Moment z';
            elseif strcmp(results_on_fig,'displacements')
                plots(1).title = 'Displacement x';
                plots(2).title = 'Displacement y';
                plots(3).title = 'Displacement z';
                plots(4).title = 'Rotation x';
                plots(5).title = 'Rotation y';
                plots(6).title = 'Rotation z';
            else
                error('Error: not recognized results_on_fig!');
            end
    
            if max(elem_id) > results.number_elements || min(elem_id) < 0
    
                error('Error: check element_id!')
    
            end
            
            figure('color',[1,1,1])
            nid = numel(elem_id);
            for i=1:nid
    
                dof = logical(results.element(elem_id(i)).node_1.node_dof);
                ttl = plots(dof).title;
                P = elems(i).Pdir;
                x1 = 0;
                x2 = results.element(i).L;
    
                subplot(nid,1,i);
                %plot([x1 x2],[P(1) P(2)],'k','linewidth',2)
                a = area([x1 x2],[P(1) P(2)]);
                a(1).FaceColor = [0 0.5 1];
                a(1).EdgeColor = [0 0 0];
                a(1).LineWidth = 2;
                a(1).FaceAlpha = 0.4;
                ax = gca; % current axes
                ax.XTick = [x1 x2];
                xlabel('Length');
                ylabel(ttl);
                title(['Element ' num2str(elem_id(i))])
      
    
            end
        end

    else

        disp('Warning: not showing results_on_fig!');

    end

end

function show_model(model,show_forces,color_nodes,color_elements,...
    color_support,holdon,box,node_num,elem_num,alpha_col)
    
    util = MSAutils;
    str2rgb = util.str2rgb;
    analysis = model.analysis_type;
    dim = str2num(analysis(end-1:end-1));

    %fig = true;

    %str_n = num2str((1:model.number_nodes)', '$\\raisebox{.5pt}{\\textcircled{\\raisebox{-.9pt} {%d}}}$');
    %str_e = num2str((1:model.number_elements)', '$\\raisebox{.5pt}{\\textcircled{\\raisebox{-.9pt} {%d}}}$');

    fid = {};
    for i=1:model.number_nodes
        x(i) = model.node(i).x;
        y(i) = model.node(i).y;
        z(i) = model.node(i).z;
        fid(i).f = sign(model.node(i).f);
    end
   
    if ~holdon
        figure('color',[1,1,1])
        hold on
    else
        hold on
    end

    minx = box(1);
    maxx = box(2);
    miny = box(3);
    maxy = box(4);
    minz = box(5);
    maxz = box(6);
    Dx = maxx - minx;
    Dy = maxy - miny;
    Dz = maxz - minz;
    Dxyz = [Dx Dy Dz]; 
    
    for i=1:model.number_elements
    
        L = model.element(i).L;
        x1=model.element(i).node_1.x;
        y1=model.element(i).node_1.y;
        z1=model.element(i).node_1.z;
        x2=model.element(i).node_2.x;
        y2=model.element(i).node_2.y;
        z2=model.element(i).node_2.z;

       rgbVector = str2rgb(color_elements);
       rgb_alpha = [rgbVector,alpha_col];
       
        if dim==1
            
            plot([x1 x2],[0 0],'Color',rgb_alpha,'linewidth',2);

            if elem_num
                str = num2str(model.element(i).id);
                text(x1+L/2, y1+Dy/20, ['$\raisebox{.5pt}{\fbox{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex', 'FontSize',14);
            end
        elseif dim==2

            plot([x1 x2],[y1 y2],'Color',rgb_alpha,'linewidth',2);
            if elem_num
                str = num2str(model.element(i).id);
                text((x1+x2)/2, (y1+y2)/2+Dy/20, ['$\raisebox{.5pt}{\fbox{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex', 'FontSize',14);
            end
        elseif dim==3
            plot3([x1 x2],[y1 y2],[z1 z2],color_elements,'linewidth',2);
        else
            error('Error: plot dimension not recognized!');
        end
    
    end



    if show_forces
         % Add arrows
        for i=1:model.number_nodes
            
            force = model.node(i).f;
            position = [model.node(i).x, model.node(i).y, model.node(i).z];
            plot_forces(force,position,Dxyz,dim,'b');
           
        end

     end
    
    for i=1:model.number_nodes
    
        plot_node(model.node(i),dim,color_nodes);
        if node_num
            str = num2str(model.node(i).id);
            x = model.node(i).x;
            y = model.node(i).y;
            z = model.node(i).z;

            if dim==1
                text(x, y+Dy/20, ['$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex','Color','m', 'FontSize',16);
            elseif dim==2
                text(x, y+Dy/20, ['$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex','Color','m', 'FontSize',16);
            elseif dim==3
                error('Error: plot dimension not recognized!');
            else

            end
        end
    
    end
 
    if dim==1

        for i=1:numel(model.support)
            % constr=model.support(i).constraint;
            plot_support_12D(model.support(i),color_support,Dx,dim);
        end

        xlim([minx - Dx/3 maxx + Dx/3]);
        ylim([minx - Dx/3 minx + Dx/3]);

        xlabel('X');
        ylabel(' ');

    elseif dim==2
    
        for i=1:numel(model.support)
             constr=model.support(i).constraint;
             plot_support(model.support(i).node,Dx,constr,color_support);
        end
         
        xlim([minx - Dx/3 maxx + Dx/3]);
        if strcmp(model.analysis_type,'beam_2D')
            ylim([-7 7]);
        else
            
            if maxy ~=miny
                ylim([miny - Dy/3 maxy + Dy/3]);
            else
                ylim([miny - Dx/3 maxy + Dx/3]);
            end
        end
         
        xlabel('X');
        ylabel('Y');
    
    else 
        xlim([minx - Dx/3 maxx + Dx/3]);
        ylim([miny - Dy/3 maxy + Dy/3]);
        zlim([minz - Dz/3 maxz + Dz/3]);
        disp('Supports in 3D not implemented yet')
    end

    axis equal;
    axis off;

end

function plot_node(node,dim,color_nodes)
    
    hinge_const = 0;
    if node.hinge
        hinge_const = 1;
    end

    if dim==1
        sc=scatter(node.x,0,'MarkerEdgeColor',color_nodes,...
                  'MarkerFaceColor',[1 1 1]*hinge_const + (1-hinge_const)*color_nodes, 'LineWidth',1.5);

    elseif dim==2

        sc=scatter(node.x,node.y,'MarkerEdgeColor',color_nodes,...
              'MarkerFaceColor',[1 1 1]*hinge_const + (1-hinge_const)*color_nodes, 'LineWidth',1.5);

    else

        sc=scatter(node.x,node.y,node.z,'MarkerEdgeColor',color_nodes,...
              'MarkerFaceColor',[1 1 1]*hinge_const + (1-hinge_const)*color_nodes, 'LineWidth',1.5);

    end

    sc.SizeData = 100;
    hold on

end

function plot_forces(force,position,Dxyz,dim,color)

    util = MSAutils;
    circular_arrow = util.circular_arrow;

    Dx = Dxyz(1);
    Dy = Dxyz(2);
    Dz = Dxyz(3);

    x = position(1);
    y = position(2);
    z = position(3);

    force(isnan(force))=0;
    if dim==1
        %idx = 1;
        F(1).force = force(1);
        F(1).type = 'F';
        num_force = 1;
    elseif dim==2
        %idx = [1 2 6];
        F(1).force = [force(1) 0];
        F(2).force = [0 force(2)];
        F(3).force = force(6);
        F(1).type = 'F';
        F(2).type = 'F';
        F(3).type = 'M';
        num_force = 3;

    elseif dim==3
        %idx = 1:6;
        num_force = 6;
        for i=1:num_force
            f0 = zeros(1,num_force);
            f0(i) = force(i);
            F(i).force = f0;
        end

        F(1).type = 'F';
        F(2).type = 'F';
        F(3).type = 'F';
        F(4).type = 'M';
        F(5).type = 'M';
        F(6).type = 'M';

    else

        error('Error: dimension cannot be larger than 3.');
        
    end

    for i=1:num_force
        
        f = F(i).force;
        fnorm = norm(f);
        if fnorm==0
            fnorm=1;
        end
        fn = f/fnorm;
        
        if strcmp(F(i).type,'F')

            if dim==1

                p1 = [x-(Dx/5)*fn 0];
                p2 = [x-Dx/50*fn  0];
                
                dp = p2-p1;    
                quiver(p1(1),p1(2),dp(1),dp(2),0,color,'LineWidth',4,'MaxHeadSize',1);  
           
                if fn ~= 0
                    text(p1(1)+fn(1)*Dx/50,p1(2)+Dx/50, sprintf(num2str(sum(abs(f))),p1),'Color',color,'FontSize',14);
                end

            elseif dim==2
                p2 = [x-(Dx/50)*fn(1)  y-(Dy/50)*fn(2)];
                p1 = [x-(Dx/5)*fn(1)  y-(Dy/5)*fn(2)];
    
                dp = p2-p1;    
                quiver(p1(1),p1(2),dp(1),dp(2),0,color,'LineWidth',2,'MaxHeadSize',1); 

                if norm(f) ~= 0
                    text(p1(1)+fn(1)*Dx/50,p1(2)+fn(1)*Dx/50, sprintf(num2str(sum(abs(f))),p1),'Color',color,'FontSize',14);
                end
            else

            end

        elseif strcmp(F(i).type,'M')
            
            if f<0
                sgn=-1;
            else
                sgn=1;
            end
            
            if f ~= 0
                radius = Dx/15; % Height from top to bottom
                centre = [x y];
                arrow_angle = 0; % Desired orientation angle in degrees
                angle = 250; % Anglebetween start and end of arrow
                direction = -sgn; % for CW enter 1, for CCW enter -1
                %colour = 'b'; % Colour of arrow
                
                head_size = 10; % Arrow head size
            
                %figHandle = figure(1); % Needs a handle
                
                hold on;    
                circular_arrow([], radius, centre, arrow_angle, angle, direction, color, head_size);
                drawnow
                
                p1 = [x+Dx/12 y+Dy/12];
                text(p1(1),p1(2), sprintf(num2str(sum(abs(f))),p1),'Color',color,'FontSize',14)
            end 
        else

        end

    end

end

function plot_support_12D(support,color_support,Dx,dim)

    hold on
    draw_fixed_support_12D(Dx,support,color_support);

end

function draw_fixed_support_12D(Dx,support,color_support)

    side_len = Dx/10;
    h0 = side_len*sin(pi/3);

    t = support.angle;
    R=[cosd(t) -sind(t);sind(t) cosd(t)];

    x0 = [support.node.x-(h0/1.5), support.node.x+(h0/1.5)];
    y0 = [0 0];

    p10 = [x0(1);y0(1)] - [mean(x0);mean(y0)];
    p20 = [x0(2);y0(2)] - [mean(x0);mean(y0)];
    p1 = R*p10+[mean(x0);mean(y0)];
    p2 = R*p20+[mean(x0);mean(y0)];

    plot([p1(1) p2(1)],[p1(2) p2(2)],color_support,'LineWidth',1);
    hold on

    %points(:,1) = p1;
    %points(:,1) = p2;

    dx = linspace(x0(1),x0(2),5);

   
    
    for i=1:numel(dx)
        p1a = [dx(i);y0(1)];
        p2a = [dx(i)+ Dx/60;y0(1)-Dx/30]; 

        m1 = [mean(dx);0];
       
        p1b = p1a - m1;
        p2b = p2a - m1;
     
        
        p11 = R*p1b + m1;
        p21 = R*p2b + m1;

        plot([p11(1) p21(1)],[p11(2) p21(2)],color_support,'LineWidth',1)
        hold on
    end

end

function plot_support(node,Dx,constr,color_support)

    side_len = Dx/15;
    h = side_len*sin(pi/3)/1.3;
    h0 = side_len*sin(pi/3);
    p = nsidedpoly(3, 'Center', [node.x ,node.y-h], 'SideLength', side_len);

    Fcol = color_support;
   

    if isequal(constr,[1 0 1 1 1 0])

        hold on
        plot([node.y-0.5*h node.y-0.5*h],[node.x-(h0/1.5), node.x+(h0/1.5)],color_support,'LineWidth',2)
        plot([node.y-1*h node.y-1*h],[node.x-(h0/1.5), node.x+(h0/1.5)],color_support,'LineWidth',2)

        sc=scatter(node.x, node.y,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1], 'LineWidth',1.5);
    
        sc.SizeData = 100;

    elseif isequal(constr,[0 1 1 1 1 0])
        
        hold on
        pg = plot(p);
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-1.8*h node.y-1.8*h],color_support,'LineWidth',2)
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 1 1 1 1 0])

        pg = plot(p);
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 1 1 1 1 1])

        hold on
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-0.5*h node.y-0.5*h],color_support,'LineWidth',2)
        pg.FaceColor = Fcol;

    elseif isequal(constr,[0 1 1 1 1 1])

        hold on
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-0.5*h node.y-0.5*h],color_support,'LineWidth',2)
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-1*h node.y-1*h],color_support,'LineWidth',2)
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 0 1 1 1 1])

        hold on
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-0.5*h node.y-0.5*h],color_support,'LineWidth',2)
        plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-1*h node.y-1*h],color_support,'LineWidth',2)
        pg.FaceColor = Fcol;

    else

        error('Error: incosistent support!')

    end
    

    
end

