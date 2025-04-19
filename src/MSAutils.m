function util = MSAutils
  util.header = @header;
  util.test_line_size = @test_line_size;
  util.pchipd = @pchipd;
  util.circular_arrow = @circular_arrow;
  util.str2rgb = @str2rgb;
  util.sort_back = @sort_back;
  util.coordinates_deformed = @coordinates_deformed;
  util.show_save_results = @show_save_results;
  util.get_force_internal = @get_force_internal;
  util.interpcolor = @interpcolor;
  util.get_flag_data = @get_flag_data;
  util.check_analysis = @check_analysis;
  util.check_variables = @check_variables;
  util.get_B = @get_B;
  util.msa2sm = @msa2sm;
  util.get_direction = @get_direction;
end

function header(fileID,input_file,analysis_type,output_file)
    fprintf(fileID,'%s\n','------------------------------------------------------------');
    fprintf(fileID,'%s\n','*            UNIVERSITY OF MINNESOTA TWIN CITIES           *');
    fprintf(fileID,'%s\n','* Department of Civil, Environmental, and Geo- Engineering *');
    fprintf(fileID,'%s\n','------------------------------------------------------------');
    fprintf(fileID,'%s\n','MSAtool: a Computational Tool for Matrix Structural Analysis');
    fprintf(fileID,'%s\n','  ');
    fprintf(fileID,'%s\n','- Purpose: education only.');
    fprintf(fileID,'%s\n','- Author: Prof. Ketson R. M. dos Santos, PhD');
    fprintf(fileID,'%s\n','- Email: dossantk@umn.edu');
    fprintf(fileID,'%s\n','- Version: 2.0.0');
    fprintf(fileID,'%s\n','  ');
    %fprintf(fileID,'%s\n','* Running: read_data');
    fprintf(fileID,'%s\n','   ');
    fprintf(fileID,'%s\n',['* Input file: ' input_file]);
    if exist(analysis_type,'var')
        fprintf(fileID,'%s\n',['* Analysis type: ' analysis_type]);
    end

    if exist(output_file,'var')
        disp(['* Output file: [' 8 [output_file ']'] 8]);
    end

end

function [opt, opt_type, coordinate_system] = get_direction(model)

    string = model.output_config.view.on_model{2};
    analysis = model.analysis_type;
  
    if strcmp(string(1),'U') || strcmp(string(1),'D')
        opt = 'displacement';
    elseif strcmp(string(1),'F') || strcmp(string(1),'P')
        opt = 'force';
    else
        error('MSAtool error: not recognized option of visualization.');
    end

    if strcmp(string,'U_x') || strcmp(string,'F_x')
        
        opt_type = 1;
        coordinate_system = 'global';

    elseif strcmp(string,'U_y') || strcmp(string,'F_y')

        opt_type = 2;
        coordinate_system = 'global';

    elseif strcmp(string,'U_xy') || strcmp(string,'F_xy')

        opt_type = 3;
        coordinate_system = 'global';

    elseif strcmp(string,'D_x') || strcmp(string,'P_x')

        opt_type = 1;
        coordinate_system = 'local';

    elseif strcmp(string,'D_y') || strcmp(string,'P_y')

        opt_type = 2;
        coordinate_system = 'local';

    elseif strcmp(string,'D_xy') || strcmp(string,'P_xy')

        opt_type = 3;
        coordinate_system = 'global';

    else
        error('MSAtool error: not recognized option of visualization.');
    end

%     if strcmp(analysis,'beam')
%     
%         if opt_type == 2
%             opt_type = 1;
%         elseif opt_type == 3
%             opt_type = 2;
%         end
% 
%     end




end

function Ps = msa2sm(P,type)

    if strcmp(type,'x')

        Ps(1) = -P(1);
        Ps(2) = P(2);

    elseif strcmp(type,'y')

        Ps(1) = P(1);
        Ps(2) = -P(2);

    elseif strcmp(type,'xy')
        Ps(1) = -P(1);
        Ps(2) = P(2);
    end

end

function B = get_B(p1,p2)

    %rng(123)
    S1 = p2-p1;
    s1 = S1/norm(S1);

    v = [s1(1);0;0];
    
  


%     t=pi/3;
% 
%     if isequal(s1,[0;1;0])
%         s2 = [-1;0;0];
%     else
%         rot = [cos(t) 0 sin(t);0 1 0;-sin(t)  0 cos(t)];
%         R = rot*s1;
%         r = R/norm(R);
%         
%         S2 = cross(s1,r);
%         s2 = S2/norm(S2);
%     end
% 
%     S3 = cross(s1,s2);
%     s3 = S3/norm(S3);
% 
% 
%     L = [s1 s2 s3];
%     
%     B = zeros(12,12);
%     
%     B(1:3,1:3) = L;
%     B(4:6,4:6) = L;
%     B(7:9,7:9) = L;
%     B(10:12,10:12) = L;

   
end

function check_variables(variable,type)

    if strcmp(type,'int')

        if floor(variable)~=variable
            error('MSA error: non integer input variable!')
        end

        if isnan(variable) || isinf(variable)
            error('MSA error: not numeric input variable!')
        end

    elseif strcmp(type,'str')

        if ~isstr(variable)
            error('MSA error: non string input variable!')
        end

    elseif strcmp(type,'double')

        if ~isnumeric(variable)
            error('MSA error: not numeric input variable!')
        end

        if isnan(variable) || isinf(variable)
            error('MSA error: not numeric input variable!')
        end

    end

end

function [dof_node, ndof_node, ndof_element]=check_analysis(analysis_type)
   


    if strcmp(analysis_type,'bar_1D')

        %dof_node = [1 0 0 0 0 0];
        dof_node = [1 0 0];
        
    elseif strcmp(analysis_type,'beam_2D')

        %dof_node = [0 1 0 0 0 2];
        dof_node = [0 1 2];

    elseif strcmp(analysis_type,'truss_2D')

        %dof_node = [1 2 0 0 0 0];
        dof_node = [1 2 0];

    elseif strcmp(analysis_type,'frame_2D')

        %dof_node = [1 2 0 0 0 3];
        dof_node = [1 2 3];

%     elseif strcmp(analysis_type,'frame_3D')
% 
%         dof_node = [1 2 3 4 5 6];
% 
%     elseif strcmp(analysis_type,'truss_3D')
% 
%         dof_node = [1 2 3 0 0 0];
% 
%     elseif strcmp(analysis_type,'grid_3D')
% 
%         dof_node = [0 2 0 1 0 3];

    else
        error('MSAtool error: analysis not not recognized!');
    end

    ndof_node = numel(dof_node(dof_node~=0));
    ndof_element = 2*ndof_node;
    
end

function line_data = get_flag_data(flag, lines_input, max_line_data)

    index = find(strcmp(flag, lines_input));
    num_lines = numel(lines_input);

    line_data = [];
    count = 1;
    bool = true;
    i=index+1;

    if numel(index)~=1
        error('MSAtool error: multiple flags with the same name found in the input file!');
    end

    while i <= num_lines && bool && count <= max_line_data
        line = lines_input{i};

        if strcmp(line(1),'*')
            bool=false;
        elseif strcmp(line(1),'%')
            % Do nothing this is a comment!  
        else
            line_data(count).cell = split(line);
            count = count + 1;
        end
        i = i + 1;
    end

end


function show_save_results(results,model,file_id)

    save_output = model.output_config.save_output;
    show_output = model.output_config.show_output;

    boolv = [save_output show_output];
    fid = [file_id 1];

    F = results.F;
    U = results.U;

    cval=1;
    for bool=boolv

        fileID = fid(cval);

        if bool
            
            fileID_ = fileID;
            if cval==2
                fileID_ = 2;
            end

            fprintf(fileID_,'%s\n','------------------------------------------------------------');
            fprintf(fileID_,'%s\n','* node_id, dof_global, nodal displacements');
    
            for i=1:results.number_nodes
                node_id = results.node(i).id;
                dof_global = results.node(i).dof_global;
                
                ndof = numel(dof_global);
                for j=1:ndof
                    fprintf(fileID,'%d %d %d \n', node_id, dof_global(j), U(dof_global(j)));
                end
            end
    
            fprintf(fileID_,'%s\n','------------------------------------------------------------');
            fprintf(fileID_,'%s\n','* node_id, dof_global, nodal forces (global)');
   
    
            for i=1:results.number_nodes
                node_id = results.node(i).id;
                dof_global = results.node(i).dof_global;
    
                ndof = numel(dof_global);
                for j=1:ndof    
                    fprintf(fileID,'%d %d %d \n', node_id, dof_global(j), F(dof_global(j)));
                end
            end
    
    
            fprintf(fileID_,'%s\n','------------------------------------------------------------');
            fprintf(fileID_,'%s\n','* node_id, dof_global, support reactions');
    
            %dof_global_constraint
            for i=1:model.number_supports
                R = results.support(i).R;
                node_id = results.support(i).node_id;
                dof_global = results.support(i).dof_global_constraint;
                
                ndof = numel(dof_global);
                
                for j=1:ndof
                    fprintf(fileID,'%d %d %d ',node_id,dof_global(j), R(j));
                    fprintf(fileID,'\n');
                end
        
            end
    
            fprintf(fileID_,'%s\n','------------------------------------------------------------');
            fprintf(fileID_,'%s\n','* Internal Forces');
            for i=1:model.number_elements
               
                P = results.element(i).P;
                fprintf(fileID,'%d %d ',i, P);
                fprintf(fileID,'\n');
        
            end
    
        end
        cval = cval + 1;
    end

end

function aux=match_data_id(data,text)

    num_data = numel(data);

    aux = data;
    for i=1:num_data
        
        nid(i) = data(i).id;

    end

    % check if ids are consistent
    if ~isequal(sort(nid),1:num_data)
        error(['Error: ' text ' IDs are not consistent!'])
    else

        aux = data(nid);
    
    end

end

function Pm = get_force_internal(element)
    
    P = element.P;
    Fnum = numel(P);

    F1 = P(1:Fnum/2);
    F2 = P(Fnum/2+1:end);

    dof_1 = element.node_1.node_dof;
    dof_2 = element.node_2.node_dof;

    sign_1 = [-1, 1, 1,-1, -1, -1];
    sign_2 = [1, -1, -1, 1, 1, 1];

    dof_sgn_1 = sign_1(dof_1==1);
    dof_sgn_2 = sign_2(dof_2==1);
    
    P1 = F1.*dof_sgn_1';
    P2 = F2.*dof_sgn_2';

    Pm = [P1;P2];

end

function yint=interpcolor(x,y,xint)
    
    
    yint = [0 0 0];
    for i=1:3
        yint(i)=interp1(x,y(:,i),xint,'linear');
    end

end

function [x,y]=coordinates_deformed(element, model_element, scale)

    x10 = model_element.node_1.x;
    y10 = model_element.node_1.y;
    z10 = model_element.node_1.z;
    x20 = model_element.node_2.x;
    y20 = model_element.node_2.y;
    z20 = model_element.node_2.z;

    x1 = element.node_1.x;
    y1 = element.node_1.y;
    z1 = element.node_1.z;
    x2 = element.node_2.x;
    y2 = element.node_2.y;
    z2 = element.node_2.z;

    dx1 = x1 - x10; 
    dy1 = y1 - y10; 
    dz1 = z1 - z10; 
    dx2 = x2 - x20; 
    dy2 = y2 - y20; 
    dz2 = z2 - z20; 

    ndof=sum(element.node_1.node_dof);
    
    Pm = get_force_internal(element);
    Ma = Pm(ndof);
    Mb = Pm(2*ndof);

    U = element.U;

    x = linspace(x1,x2,100);
    if element.node_1.node_dof(6)==0
        y = linspace(y1,y2,100);
    else
        
        %scale=0.2/Tmax;
        ndof = numel(element.dof_global);
        t1 = U(ndof/2);
        t2 = U(ndof);

        if x1==x2 && t1==0 && t2==0
            y = linspace(y1,y2,100);

        else
            L = element.L;
            theta = element.angle;
            R = [cosd(theta) -sind(theta);
                 sind(theta) cosd(theta)];
           
            if abs(theta)~=90
                Ta = scale*[tan(t1) tan(t2)] + [tand(theta) tand(theta)];
                x0 = linspace(x1, x2, 100);
                yp = pchipd([x1 x2],[y1 y2],Ta, x0);
                zr = [x0; yp];
       
            else
                Ta = scale*[tan(t1) tan(t2)];
                x0 = linspace(0, L, 100);
                yp = pchipd([0 L],[0 0],Ta,x0);
    
                c1 = ([dx2;dy2] - [dx1;dy1]).*linspace(0,1,100) + [dx1;dy1];
                %c1 = ([dx2;0] - [dx1;0]).*linspace(0,1,100) + [dx1;0];
                z = [x0;yp];
                zr = R*z + c1 + [x10;y10];
            end
 
            x = zr(1,:);
            y = zr(2,:);
 

        end
    end

end

function test_line_size(line_data, size, text)

    num_line = sum(~cellfun('isempty',line_data));
    
    if num_line ~= size
    
        error(['Error: Input data for '...
            text ': ' num2str(size)...
            ' numbers/data are expected per line, but only '...
            num2str(num_line) ' are provided in at least one of them.' ]);

    end

end

function circular_arrow(figHandle, radius, centre, arrow_angle, angle, direction, colour, head_size, head_style)
    % This is a function designed to draw a circular arrow onto the current
    % figure. It is required that "hold on" must be called before calling this
    % function. 
    %
    % The correct calling syntax is:
    %   circular_arrow(height, centre, angle, direction, colour, head_size)
    %   where:
    %       figHandle - the handle of the figure to be drawn on.
    %       radius - the radius of the arrow. 
    %       centre - a vector containing the desired centre of the circular
    %                   arrow.
    %       arrow_angle - the desired orientation angle of the circular arrow.
    %                   This is measured in degrees counter-clockwise 
    %       angle - the angle between starting and end point of the arrow in
    %                   degrees.
    %       direction - variable set to determine format of arrow head. Use 1
    %                   to get a clockwise arrow, -1 to get a counter clockwise
    %                   arrow, 2 to get a double headed arrow and 0 to get just
    %                   an arc. 
    %       colour (optional) - the desired colour of the arrow, using Matlab's
    %                   <a href="matlab:
    %                   web('https://au.mathworks.com/help/matlab/ref/colorspec.html')">Color Specification</a>. 
    %       head_size (optional) - the size of the arrow head.
    %       head_style (optional) - the style of the arrow head.
    %                   For more information, see <a href="matlab: 
    %                   web('http://au.mathworks.com/help/matlab/ref/annotationarrow-properties.html#property_HeadStyle')">Annotation Arrow Properties</a>.
    %Ensure proper number of arguments

  
    if (nargin < 6)||(nargin > 9)
        error(['Wrong number of parameters '...
            'Enter "help circular_arrow" for more information']);
    end
    % arguments 7, 8 and 9 are optional,
    if nargin < 9
       head_style = 'vback2';
    end
    if nargin < 8
       head_size = 10;
    end
    if nargin < 7
       colour = 'k';
    end
    % display a warning if the headstyle has been specified, but direction has
    % been set to no heads
    if nargin == 9 && direction == 0
        warning(['Head style specified, but direction set to 0! '...
            'This will result in no arrow head being displayed.']);
    end
        
    % Check centre is vector with two points
    [m,n] = size(centre);
    if m*n ~= 2
        error('Centre must be a two element vector');
    end
    arrow_angle = deg2rad(arrow_angle); % Convert angle to rad
    angle = deg2rad(angle); % Convert angle to rad
    xc = centre(1);
    yc = centre(2);
    % Creating (x, y) values that are in the positive direction along the x
    % axis and the same height as the centre
    x_temp = centre(1) + radius;
    y_temp = centre(2);
    % Creating x & y values for the start and end points of arc
    x1 = (x_temp-xc)*cos(arrow_angle+angle/2) - ...
            (y_temp-yc)*sin(arrow_angle+angle/2) + xc;
    x2 = (x_temp-xc)*cos(arrow_angle-angle/2) - ...
            (y_temp-yc)*sin(arrow_angle-angle/2) + xc;
    x0 = (x_temp-xc)*cos(arrow_angle) - ...
            (y_temp-yc)*sin(arrow_angle) + xc;
    y1 = (x_temp-xc)*sin(arrow_angle+angle/2) + ...
            (y_temp-yc)*cos(arrow_angle+angle/2) + yc;
    y2 = (x_temp-xc)*sin(arrow_angle-angle/2) + ... 
            (y_temp-yc)*cos(arrow_angle-angle/2) + yc;
    y0 = (x_temp-xc)*sin(arrow_angle) + ... 
            (y_temp-yc)*cos(arrow_angle) + yc;
    % Plotting twice to get angles greater than 180
    i = 1;
    % Creating points
    P1 = struct([]);
    P2 = struct([]);
    P1{1} = [x1;y1]; % Point 1 - 1
    P1{2} = [x2;y2]; % Point 1 - 2
    P2{1} = [x0;y0]; % Point 2 - 1
    P2{2} = [x0;y0]; % Point 2 - 1
    centre = [xc;yc]; % guarenteeing centre is the right dimension
    n = 1000; % The number of points in the arc
    v = struct([]);
       
    while i < 3
        v1 = P1{i}-centre;
        v2 = P2{i}-centre;
        c = det([v1,v2]); % "cross product" of v1 and v2
        a = linspace(0,atan2(abs(c),dot(v1,v2)),n); % Angle range
        v3 = [0,-c;c,0]*v1; % v3 lies in plane of v1 and v2 and is orthog. to v1
        v{i} = v1*cos(a)+((norm(v1)/norm(v3))*v3)*sin(a); % Arc, center at (0,0)
        plot(v{i}(1,:)+xc,v{i}(2,:)+yc,'Color', colour,'linewidth',1); % Plot arc, centered at P0
        
        i = i + 1;
    end
  
    position = struct([]);
    % Setting x and y for CW and CCW arrows
    if direction == 1
        position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
    elseif direction == -1
        position{1} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];
    elseif direction == 2
        position{1} = [x2 y2 x2-(v{2}(1,2)+xc) y2-(v{2}(2,2)+yc)];
        position{2} = [x1 y1 x1-(v{1}(1,2)+xc) y1-(v{1}(2,2)+yc)];  
    elseif direction == 0
        % Do nothing
    else
        error('direction flag not 1, -1, 2 or 0.');
    end
    % Loop for each arrow head
    i = 1;
    while i < abs(direction) + 1
        h=annotation('arrow'); % arrow head
        set(h,'parent', gca, 'position', position{i}, ...
            'HeadLength', head_size, 'HeadWidth', head_size,...
            'HeadStyle', head_style, 'linestyle','none','Color', colour);
        i = i + 1;
    end

end

function rgbVector = str2rgb(colorString)
%STR2RGB   Converts a string representation of a color to an RGB triple.
%   X = STR2RGB(STR) converts the string STR into a three-element row
%   vector X (an RGB triple). STR can be one of three string
%   representations of colors that MATLAB uses (see ColorSpec help): a full
%   color name ('yellow'), a single character ('y'), or a string of numbers
%   within the range of 0 and 1 ('[1 1 0]' or '1 1 0').
%
%   If the string STR does not represent a valid string representation of a
%   color, STR2RGB(STR) returns NaN.
%
%   NOTE: STR2RGB does not use eval.
%
%   Examples:
%      str2rgb('yellow')        returns [1 1 0]
%      str2rgb('y')             returns [1 1 0]
%      str2rgb('[1 1 0]')       returns [1 1 0]
%      str2rgb('1 1 0')         returns [1 1 0]
%      str2rgb('[1; 1; 0]')     returns [1 1 0]
%      str2rgb('[0 0.5 0.91]')  returns [0 0.5000 0.9100]
%      str2rgb('purple')        returns NaN
%      str2rgb('[1 2]')         returns NaN
% Author: Ken Eaton
% Last modified: 4/2/08
%--------------------------------------------------------------------------
  if (nargin == 0),
    error('str2rgb:notEnoughInputs','Not enough input arguments.');
  end
  if (~ischar(colorString)),
    error('str2rgb:badArgumentType',...
          'Input argument should be of type char.');
  end
  expression = {'^\s*yellow\s*$','^\s*magenta\s*$','^\s*cyan\s*$',...
                '^\s*red\s*$','^\s*green\s*$','^\s*blue\s*$',...
                '^\s*white\s*$','^\s*black\s*$','^\s*y\s*$','^\s*m\s*$',...
                '^\s*c\s*$','^\s*r\s*$','^\s*g\s*$','^\s*b\s*$',...
                '^\s*w\s*$','^\s*k\s*$','[\[\]\;\,]'};
  replace = {'[1 1 0]','[1 0 1]','[0 1 1]','[1 0 0]','[0 1 0]',...
             '[0 0 1]','[1 1 1]','[0 0 0]','[1 1 0]','[1 0 1]',...
             '[0 1 1]','[1 0 0]','[0 1 0]','[0 0 1]','[1 1 1]',...
             '[0 0 0]',' '};
  rgbVector = sscanf(regexprep(colorString,expression,replace),'%f').';
  if ((numel(rgbVector) ~= 3) || any((rgbVector < 0) | (rgbVector > 1))),
    rgbVector = nan;
  end
  
end

function pp = pchipd(x,y,d,xx)
    %PCHIPD  Piecewise Cubic Hermite Interpolating Polynomial with Derivatives.
    %   PP = PCHIPD(X,Y,D) provides the piecewise cubic polynomial which
    %   interpolates values Y and derivatives D at the sites X.  This is meant
    %   to augment the built-in Matlab function PCHIP, which does not allow the
    %   user to specify derivatives.
    %  
    %   X must be a vector.
    %
    %   If Y and D are vectors, then Y(i) and D(i) are the value and derivative
    %   to be matched at X(i).
    %
    %   If Y and D are matrices, then size(Y,2) == size(D,2) == length(X).
    %   Also, size(Y,1) == size(D,1).  Use this for interpolating vector valued
    %   functions.
    %
    %   YY = PCHIPD(X,Y,D,XX) is the same as YY = PPVAL(PCHIPD(X,Y,D),XX), thus
    %   providing, in YY, the values of the interpolant at XX.
    %
    %   Example comparing SPLINE, PCHIP, and PCHIPD
    %     a = -10;
    %     b = 10;
    %     x = linspace(a,b,7); 
    %     f = @(x) 1./(1+exp(-x));  % logistic function
    %     df = @(x) f(x).*(1-f(x)); % derivative of the logistic function
    %     t = linspace(a,b,50);
    %     r = f(t);
    %     p = pchip(x,f(x),t);
    %     s = spline(x,f(x),t);
    %     q = pchipd(x,f(x),df(x),t);
    %     plot(t,r,'k',x,f(x),'o',t,p,'-',t,s,'-.',t,q,'--')
    %     legend('true','data','pchip','spline','pchipd',4)
    %
    %   See also INTERP1, SPLINE, PCHIP, PPVAL, MKPP, UNMKPP.
    %
    %
    % 2010-10-04 (nwh) first version
    %
    % check inputs
    % x must be a vector
    if ~isvector(x)
      error('pchipd:input_error','x must be a vector of length > 2.')
    end
    % get size and orient
    n = length(x);
    x = x(:);
    % make sure x is long enough, we can't construct an interpolating
    % polynomial with just one point
    if n < 2
      error('pchipd:input_error','x must be a vector of length > 2.')
    end
    % check y and d
    if isvector(y) && isvector(d) && length(y) == n && length(d) == n
      % orient
      y = y(:);
      d = d(:);
      m = 1;
    elseif size(y,2) == n && size(d,2) == n && size(y,1) == size(d,1)
      m = size(y,1);
      y = y';
      d = d';
    else
      error('pchipd:input_error','y and d must be vectors or matrices of same size with length(x) columns.')
    end
    % sort breaks & data if needed
    if ~issorted(x)
      [x x_ix] = sort(x);
      if m == 1
        y = y(x_ix);
        d = d(x_ix);
      else
        y = y(x_ix,:);
        d = d(x_ix,:);
      end
    end
    % compute coefficients
    coef = zeros(m,n-1,4);
    dx = diff(x);
    for i = 1:m
      dy = diff(y(:,i));
      coef(i,:,4) = y(1:end-1,i)';
      coef(i,:,3) = d(1:end-1,i)';
      coef(i,:,2) = 3*dy./(dx.^2) - (2*d(1:end-1,i)+d(2:end,i))./dx;
      coef(i,:,1) = -2*dy./(dx.^3) + (d(1:end-1,i)+d(2:end,i))./(dx.^2);
    end
    % create the piecewise polynomial structure
    pp = mkpp(x,coef,m);
    % if user requests evaluations
    if nargin > 3
      pp = ppval(pp,xx);
    end
    
end

function out = sort_back( data, ind, dim )
% SORT_BACK sort back data to original order
% ind is the indexes obtained from sorting
% dim is the sorted dimension of the data (assumed to be 1 if not specified)
% Ex:
% y = randn(3,4,2);
% [y,ind] = sort(y,2);
% do stuff with sorted y...
% y2 = sort_back( y, ind, 2 );
% 
% Works on arrays of any dimension
% Also works on cellstrings (vectors)
% 
% C = {'hello' 'yes' 'no' 'goodbye'};
% [C,ind] = sort(C);
% C2 = sort_back(C,ind);
%
% See also SORT
%Author Ivar Eskerud Smith
if size(ind)~=size(data)
    error('Different size of indexes and input data');
end
if iscell(data)
    if ~any(size(data)==1)
        error('Only vectors are supported in cell sorting/back-sorting');
    end
    out=cell(size(data));
    out(ind) = data;
    return;
end
if ~isnumeric(data) || ~isnumeric(ind)
    error('Inputs have to be numeric or cell');
end
n=ndims(ind);
if ~exist('dim','var')
    dim=1;
end
if dim>n
    error('Specified sorted dimension must be within array bounds');
end
%shift array so that the sorted dim is the first dimension
if dim~=1
    sortInd=1:1:n;sortInd(1)=dim;sortInd(dim)=1;
    data = permute(data,sortInd);
    ind = permute(ind,sortInd);
end
inds = repmat({1},1,n);inds{1}=':';
%if ~issorted( data(inds{:}) )
%    warning('The input data is not sorted along the specified dimension');
%end
s = size(ind);
nData = numel(data);
inds = repmat({1},1,n);
inds(1:2)={':',':'};
shiftSize = s(1)*s(2);
out=nan(size(data));
%loop all 2d arrays within nd-array
for k=1:prod(s(3:end))
    tmpdata = data(inds{:});
    tmpind  = ind(inds{:});
    %data is shifted so that the sorted dim = 1
    for i=1:numel(tmpdata(1,:))
        out(tmpind(:,i),i) = tmpdata(:,i);
    end
    if n>2
        %shift to next 2d array within nd-array
        shiftInds = mod((1:nData)-shiftSize-1,nData)+1;
        out=reshape(out(shiftInds),s);
        data=reshape(data(shiftInds),s);
        ind=reshape(ind(shiftInds),s);
    end
end
%permute back to original order
sortInd=1:1:ndims(out);sortInd(1)=dim;sortInd(dim)=1;
out = permute(out,sortInd);

end
