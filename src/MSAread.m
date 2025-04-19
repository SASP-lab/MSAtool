function data = MSAread
    data.msa_read_data=@msa_read_data;
end

%% ========================== MSA READ DATA ===============================

function model=msa_read_data(varargin)

    util = MSAutils;

    % Set functions.
    %test_line_size = util.test_line_size;
    get_flag_data = util.get_flag_data;
    %check_analysis = util.check_analysis;
    check_variables = util.check_variables;
    model = {};
    
    % Get boolean variables and the name of the input file.
    input_file=varargin{1};
    verbose=varargin{2};
    homedir=varargin{3};
    inputdir=varargin{4};
    %MSAtooldir=varargin{5};
    
    % Store the name of the input file in model.
    model.input_file = input_file;
    input_path = fullfile(inputdir,input_file);
    model.input_path = input_path;

    % Check if input file exists.
    if exist(input_path,'file')

        % Print the data on command window.
        if verbose
            util.header(1,model.input_file,[],[]);
        end
        file = input_path;
        fid = fopen(file,'r');
    else
         % File does not exist.
         cd(homedir)
         error('MSA tool error: input file does not exist!');
    end
    
    % Get the lines of input_file.
    lines_scan = textscan(fid,'%s','delimiter','\r');
    lines_input = lines_scan{1};

    % Read the analysis type.
    line_data_0 = get_flag_data('*ANALYSIS', lines_input, 1);
    analysis_type = line_data_0(1).cell{1};

    line_data_1 = get_flag_data('*NUMBER_MATERIALS_SECTIONS', lines_input, 1);
    number_materials = str2double(line_data_1(1).cell{1});
    number_sections = str2double(line_data_1(1).cell{2});

    check_variables(number_materials,'int');
    check_variables(number_sections,'int');
    if number_materials < 1 || number_sections < 1
        error('MSAtool error: include at least one material/section.');
    end

    line_data_mat = get_flag_data('*MATERIALS', lines_input, number_materials);
    line_data_sec = get_flag_data('*SECTIONS', lines_input, number_sections);

    material = struct('id',0,'E',0);
    for i=1:number_materials
        material(i).id = str2double(line_data_mat(i).cell{1});
        material(i).E = str2double(line_data_mat(i).cell{2});
    end
    section = struct('id',0,'A',0,'I',0);
    for i=1:number_sections
        section(i).id = str2double(line_data_sec(i).cell{1});
        section(i).A = str2double(line_data_sec(i).cell{2});
        section(i).I = str2double(line_data_sec(i).cell{3});
    end

    % Read the number of nodes, elements, and supports.
    line_data_2 = get_flag_data('*NUMBER_NODES_ELEMENTS_SUPPORTS', lines_input, 1);
    number_nodes = str2double(line_data_2(1).cell{1});
    number_elements = str2double(line_data_2(1).cell{2});
    number_supports = str2double(line_data_2(1).cell{3});

    check_variables(number_nodes,'int');
    check_variables(number_elements,'int');
    check_variables(number_supports,'int');

    if number_elements > (number_nodes*(number_nodes-1)/2)
        error('MSAtool error: the number of elements is too large.');
    elseif number_elements < number_nodes-1
        error('MSAtool error: the number of elements is too small.');
    end

    if number_supports > number_nodes
        error('MSAtool error: the number of support is larger than the number of nodes.');
    end

    % Read the coordinates of each node.
    line_data_3 = get_flag_data('*NODES_COORDINATES', lines_input, number_nodes);

    % Read the forces applied in each node.
    line_data_4 = get_flag_data('*NODES_FORCES', lines_input, number_nodes);

    % Read the prescribed displacements in the nodes.
    line_data_5 = get_flag_data('*NODES_DISPLACEMENTS', lines_input, number_nodes);
    
    % Create a structure for node.
    node = struct('id',0,'x',0,'y',0,'f',zeros(3,1),'u',zeros(3,1));
    %num_dof_global = number_nodes * ndof_node;

    %F = zeros(num_dof_global,1);
    %U = nan(num_dof_global,1);

    idnodes = zeros(1,number_nodes);
    for i=1:number_nodes
        node_data = str2double(line_data_3(i).cell);
        node_id = node_data(1);

        node(node_id).id = node_data(1); 
        node(node_id).x = node_data(2);
        node(node_id).y = node_data(3); 

        idnodes(i) = node_id;
        
        node_forces = str2double(line_data_4(i).cell);
        node_displacements = str2double(line_data_5(i).cell);
        
        check_variables(node_data(1),'int');
        check_variables(node_data(2),'double');
        check_variables(node_data(3),'double');
        check_variables(node_data(4),'int');
        
        check_variables(node_forces(1),'int');
        check_variables(node_forces(2),'double');
        check_variables(node_forces(3),'double');
        check_variables(node_forces(4),'double');
        
        check_variables(node_displacements(1),'int');
        check_variables(node_displacements(2),'double');
        check_variables(node_displacements(3),'double');
        check_variables(node_displacements(4),'double');

        node_force_id = node_forces(1);
        node_displacements_id = node_displacements(1);
        force = node_forces(2:4);
        displacements = node_displacements(2:4);

        node(node_force_id).f = force;
        node(node_displacements_id).u_p =  displacements;

        if strcmp(analysis_type,'truss')
            node(node_displacements_id).hinge = 1;
        else
            node(node_displacements_id).hinge = node_data(4);
        end

    end

    if ~isequal(sort(idnodes),1:number_nodes)
        error('MSAtool: error in node id.')
    end

    % Read the elements.
    line_data_6 = get_flag_data('*ELEMENTS', lines_input, number_elements);
    line_data_6h = get_flag_data('*ELEMENTS_HINGE', lines_input, number_elements);
    line_data_6r = get_flag_data('*ELEMENTS_RIGIDITY', lines_input, number_elements);

    % Create a structure for elements (id, E, G, A, Iy, Iz, node_1, node_2).
    %K = zeros(num_dof_global,num_dof_global);
    element = struct('id',0,'material',0,'section',0,'node_1',0,'node_2',0,'hinge_1',0,'hinge_2',0,'rigidity',0, 'rigid', 0);
    idelem = zeros(1,number_elements);
    Emod = zeros(number_elements,1);
    for i=1:number_elements
        element_data = str2double(line_data_6(i).cell);
     
        check_variables(element_data(1),'int');
        check_variables(element_data(2),'int');
        check_variables(element_data(3),'int');
        check_variables(element_data(4),'int');
        check_variables(element_data(5),'int');
        check_variables(element_data(6),'int');

        element_id = element_data(1);
        idelem(i) = element_id;
    
        element(element_id).id = element_id; 
        element(element_id).material = material(element_data(2));
        element(element_id).section = section(element_data(3));
        element(element_id).rigid = logical(element_data(4));
        element(element_id).node_1 = node(element_data(5));
        element(element_id).node_2 = node(element_data(6));
       
        element_data = str2double(line_data_6h(i).cell);
        element_id_h = element_data(1);
        element(element_id_h).hinge = logical(element_data(2:3));

        element_data = str2double(line_data_6r(i).cell);
        element_id_r = element_data(1);
        element(element_id_r).rigidity = logical(element_data(2:3));
        
        Emod(element_id) = element(element_id).material.E;
        %element(element_id).axially_rigid = element_data(8);
        
        %element(element_id).number_dof = ndof_element;

        %element(element_id).dof_global = [element(element_id).node_1.dof_global...
        %                            element(element_id).node_2.dof_global];


       
    end


    maxE = max(Emod);
    for i=1:number_elements
        element_data = str2double(line_data_6(i).cell);
        element_id = element_data(1);

        if element(element_id).rigid
            element(element_id).material.E = maxE*1e10;
        end

    end

    if ~isequal(sort(idelem),1:number_elements)
        error('MSAtool: error in element id.')
    end

    % Read the elements.
    line_data_7 = get_flag_data('*SUPPORTS', lines_input, number_supports);
   
    % Create a structure for elements (id, E, G, A, Iy, Iz, node_1, node_2).
    support = struct('node_id',0,'constraint',0,'angle',0);
   
    if numel(line_data_7)~=number_supports
        error('MSAtool: number of support data is inconsistent with the defined number of supports.')
    end

    for i=1:number_supports
        support_data = str2double(line_data_7(i).cell);

        check_variables(support_data(1),'int');
        check_variables(support_data(2),'int');
        check_variables(support_data(3),'int');
        check_variables(support_data(4),'int');
        check_variables(support_data(5),'double');

        node_id = support_data(1);
        support(i).node_id = node_id; 
        support(i).constraint = support_data(2:4);
        support(i).angle = support_data(5);

    end

    line_data_8 = get_flag_data('*OUTPUT', lines_input, 1);
    show_output = str2double(line_data_8(1).cell{1});
    save_output = str2double(line_data_8(1).cell{2});
    visual_bool = str2double(line_data_8(1).cell{3}); 

    check_variables(show_output,'int');
    check_variables(save_output,'int');
    check_variables(visual_bool,'int');

    output_config = struct('show_output',false,'save_output',false,...
                           'visualization',false);
    if show_output ~= 0 && show_output ~= 1
        error('MSAtool error: show output must be either 0 or 1.')
    else
        output_config.show_output = logical(show_output);
    end
              
    if save_output ~= 0 && save_output ~= 1
        error('MSAtool error: save output must be either 0 or 1.')
    else
        output_config.save_output = logical(save_output);
    end

    if visual_bool ~= 0 && visual_bool ~= 1
        error('MSAtool error: visualization must be either 0 or 1.')
    else

        visual_bool = logical(visual_bool);
        output_config.visualization = visual_bool;

        if ~visual_bool
            fprintf(2,'Warning MSAtool: no graphical output will be produced!\n')
        else 
            
            line_data_view = get_flag_data('*VIEW_OPTIONS', lines_input, 1);
            line_data_results = get_flag_data('*VIEW_RESULTS', lines_input, 2);
            line_data_diagrams = get_flag_data('*INTERNAL_FORCES_DIAGRAMS', lines_input, 2);

            view_opt = str2double(line_data_view.cell);
      
            for jj=1:5
                check_variables(view_opt(jj),'int');
            end
      
            % show id of nodes/elements, show applied forces, show reactions, show numerical, show deformed shape
            output_config.view.id = view_opt(1);
            output_config.view.forces = view_opt(2);
            output_config.view.reactions = view_opt(3);
            output_config.view.numerical = view_opt(4);
            output_config.view.deformed = view_opt(5);

            output_config.view.on_model{1} = str2double(line_data_results(1).cell{1});
            check_variables(output_config.view.on_model{1},'int');
            output_config.view.on_model{1} = logical(str2double(line_data_results(1).cell{1}));
            output_config.view.on_model{2} = line_data_results(1).cell{2};

            output_config.view.plot{1} = str2double(line_data_results(2).cell{1});
            check_variables(output_config.view.plot{1},'int');
            output_config.view.plot{1} = logical(str2double(line_data_results(2).cell{1}));
            output_config.view.plot{2} = line_data_results(1).cell{2};
            
            for i=3:numel(line_data_results(2).cell)
                value(i-2) = str2double(line_data_results(2).cell{i});
            end
   
            output_config.view.plot{3} = value;

            output_config.diagrams.on_model{1} = str2double(line_data_diagrams(1).cell{1});
            check_variables(output_config.diagrams.on_model{1},'int');
            output_config.diagrams.on_model{1} = logical(str2double(line_data_diagrams(1).cell{1}));
            output_config.diagrams.on_model{2} = line_data_diagrams(1).cell{2};

            output_config.diagrams.plot{1} = str2double(line_data_diagrams(2).cell{1});
            check_variables(output_config.diagrams.plot{1},'int');
            output_config.diagrams.plot{1} = logical(str2double(line_data_diagrams(2).cell{1}));
            output_config.diagrams.plot{2} = line_data_diagrams(1).cell{2};
            
            for i=3:numel(line_data_diagrams(2).cell)
                value(i-2) = str2double(line_data_diagrams(2).cell{i});
            end
            output_config.diagrams.plot{3} = value;

        end
 
    end
   
    model.analysis_type = analysis_type;
    model.number_nodes = number_nodes;
    model.number_elements = number_elements;
    model.number_supports = number_supports;
    model.node = node;
    model.element = element;
    model.support = support;
    model.output_config = output_config;  
 
end



