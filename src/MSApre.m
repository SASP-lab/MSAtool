function preprocessor = MSApre
    preprocessor.msa_pre=@msa_pre;
end

%% ========================== MSA PREPROCESSOR ===============================

function model=msa_pre(model)

    %util = MSAutils;

    % Process each element to compute the corresponding stiffness matrices
    % and compute the vectors of nodal forces and displacements.

    num_nodes = model.number_nodes;
    num_elem = model.number_elements;
    num_supports = model.number_supports;
    analysis = model.analysis_type;

    if strcmp(analysis,'bar')
        dofs = 1;
    elseif strcmp(analysis,'beam')
        dofs = [2 3];
    elseif strcmp(analysis,'truss')
        dofs = [1 2];
    elseif strcmp(analysis,'frame')
        dofs = [1 2 3];
    end
    ndof = numel(dofs);
    dof_loc = 1:ndof;
    model.dof_node = dofs;

    node = model.node;
    released = zeros(num_nodes,1);
    F = zeros(ndof*num_nodes,1);
    Up = zeros(ndof*num_nodes,1);

    for i=1:num_nodes

        nid = model.node(i).id;
        node(nid).dof = ndof*(model.node(nid).id-1) + dof_loc;
        node(nid).num_hinge = 0;
        node(nid).num_elements = 0;
        node(nid).directions = dofs;
        
        fnode = node(nid).f;

        F(node(nid).dof) = fnode(dofs);
        model.node(nid).dof = node(nid).dof;
        
        u_p = zeros(ndof,1);
        Up(node(nid).dof) = u_p;

    end

 
    F_ext = F; % applied forces.
    K = zeros(num_nodes*ndof,num_nodes*ndof);
    for i=1:num_elem

        elem = model.element(i);

        if elem.node_1.hinge
            elem.hinge(1) = true;
        end

        if elem.node_2.hinge
            elem.hinge(2) = true;
        end

        id1 = elem.node_1.id;
        id2 = elem.node_2.id;

        elem.node_1 = node(id1);
        elem.node_2 = node(id2);

        dof_global = [node(id1).dof node(id2).dof];
        elem.dof_global = dof_global;

        f = [node(id1).f;node(id2).f];
        %u_p = [node(id1).u_p;node(id2).u_p];
   
        %dof_global_hinge = [node(id1).dof_hinge node(id2).dof_hinge]; 
        %logical(dof_global_hinge)
        
        ct = double(elem.hinge)';
       
        node(id1).num_hinge = node(id1).num_hinge + ct(1);
        node(id2).num_hinge = node(id2).num_hinge + ct(2);
        node(id1).num_elements = node(id1).num_elements + 1;
        node(id2).num_elements = node(id2).num_elements + 1;
        
        [Kl, B, L, f_out,Kl_full,theta] = elemental_matrices_frame(elem,f);
        
        D = zeros(2*ndof,num_nodes*ndof);

        for j=1:2*ndof
            D(j,dof_global(j))=1;
        end

        Kg = B'*Kl*B;
        Kg_full = B'*Kl_full*B;
        
        model.element(i).node_1 = node(id1);
        model.element(i).node_2 = node(id2);
        model.element(i).L = L;
        model.element(i).K_local_fix = Kl_full;
        model.element(i).K_local = Kl;
        model.element(i).K_global = Kg;
        model.element(i).K_global_fix = Kg_full;
        model.element(i).B = B;
        model.element(i).D = D;
        model.element(i).dof = dof_global;
        model.element(i).f = f_out;
        model.element(i).angle = theta;
        
        K = K + D'*Kg*D;
   
    end
 

    for i=1:num_nodes

        if strcmp(model.analysis_type,'frame')

            nid = model.node(i).id;
            if node(nid).hinge
                released(i) = node(nid).dof(3);
            end
    
            if node(nid).num_hinge==node(nid).num_elements
                released(i) = node(nid).dof(3);
            end

        elseif strcmp(model.analysis_type,'beam')

            nid = model.node(i).id;
            if node(nid).hinge
                released(i) = node(nid).dof(2);
            end
    
            if node(nid).num_hinge==node(nid).num_elements
                released(i) = node(nid).dof(2);
            end

        end

    end

    U = nan(ndof*num_nodes,1);
    for i=1:num_supports
        sup = model.support(i);
        id = sup.node_id;
        model.support(i).node = node(id);
        constraint = logical(sup.constraint(dofs));

        F(node(id).dof(constraint)) = nan(size(F(node(id).dof(constraint))));
        U(node(id).dof(constraint)) = Up(node(id).dof(constraint));

    end
  
 
%     released(released==0)=[];
%  
%     if sum(released)>0
%         K(:,released) =[];
%         K(released,:) =[];
%         F(released) = [];
%         U(released) = [];
%     end

    model.F_ext = F_ext;
    model.K = K;
    model.F = F;
    model.U = U;
    model.node = node;
    model.released = released;
    
end

function [K_frame,B, L, f_out, K_frame_full, theta] = elemental_matrices_frame(element, f)

    %util = MSAutils;

    % Material properties.
    E = element.material.E;
    A = element.section.A;
    I = element.section.I;
    
    %element.node_1.directions
    %element.dof_global
    % Nodes positions.
    posi_1(1) = element.node_1.x;
    posi_1(2) = element.node_1.y;
    posi_2(1) = element.node_2.x;
    posi_2(2) = element.node_2.y;

    % Lenght.
    L = norm(posi_2 - posi_1);

    % Angle
    dposi = posi_2 - posi_1;
    theta = atan2d(dposi(2),dposi(1));
    c = cosd(theta);
    s = sind(theta);


    B = [c s 0 0 0 0;
        -s c 0 0 0 0;
         0 0 1 0 0 0;
         0 0 0 c s 0;
         0 0 0 -s c 0;
         0 0 0 0 0 1];

    if L<=0
        error('MSAtool error: cannot compute the length of an element.');
    end

    EA=E*A;
    EI = E*I;

    Kt = (EA/L)*[1 0 0 -1 0 0;
                 0 0 0 0 0 0;
                 0 0 0 0 0 0;
                 -1 0 0 1 0 0;
                 0 0 0 0 0 0;
                 0 0 0 0 0 0];

    id_beam = [2 3 5 6];
    fb = f(id_beam);
    fm = fb;

     Kbfull = (EI/(L^3))*[0 0 0 0 0 0;
                     0 12 6*L 0 -12 6*L;
                     0 6*L 4*(L^2) 0 -6*L 2*(L^2);
                     0 0 0 0 0 0;
                     0 -12 -6*L 0 12 -6*L;
                     0 6*L 2*(L^2) 0 -6*L 4*(L^2)];

    if element.hinge(1) && ~element.hinge(2)

        Kb = (3*EI/(L^3))*[0 0 0 0 0 0;
                           0 1 0 0 -1 L;
                           0 0 0 0 0 0;
                           0 0 0 0 0 0;
                           0 -1 0 0 1 -L;
                           0 L 0 0 -L (L^2)];

        fm = fb - (1/2*L)*[3;1;-3;L]*fb(2); 
       

    elseif ~element.hinge(1) && element.hinge(2)

        Kb = (3*EI/(L^3))*[0 0 0 0 0 0;
                           0 1 L 0 -1 0;
                           0 L (L^2) 0 -L 0;
                           0 0 0 0 0 0;
                           0 -1 -L 0 1 0;
                           0 0 0 0 0 0];

        fm = fb - (1/2*L)*[3;L;-3;1]*fb(4); 

    elseif element.hinge(1) && element.hinge(2)

        Kb = (EI/(L^3))*zeros(6,6);
        ct = ((fb(2)+fb(4))/L);
        fm = [fb(1) - ct;0;fb(3) + ct;0];
        
    else

        Kb = Kbfull;
    end
    

    K_frame = Kt + Kb;
    K_frame_full = Kt + Kbfull;

   
    rig = ones(3,1).*1e15;
    %rig = nan(3,1);
    rig(~element.rigidity) = 1;

    Kr = [rig(1) 0 0 rig(1) 0 0;
          0 rig(2) rig(2) 0 rig(2) rig(2);
          0 rig(2) rig(2) 0 rig(2) rig(2);
          rig(1) 0 0 rig(1) 0 0;
          0 rig(2) rig(2) 0 rig(2) rig(2);
          0 rig(2) rig(2) 0 rig(2) rig(2)];
    
    f(id_beam) = fm;
    f_out = f;
    K_frame = K_frame.*Kr;
    K_frame_full = K_frame_full.*Kr;

    %K_frame,B, L, f_out, K_frame_full, theta
    idx = [element.node_1.directions element.node_1.directions+3];
    K_frame = K_frame(idx,idx);
    B = B(idx,idx);
    f_out = f_out(idx);
    K_frame_full = K_frame_full(idx,idx);
    
end


