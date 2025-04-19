function viewer = MSAviewer
    viewer.msa_viewer=@msa_viewer;
    viewer.msa_viewer_user=@msa_viewer_user;
end

%% ========================== MSA RESULTS =================================

%function msa_results(model,results,results_on_fig,fig_model)

function msa_viewer_user(model,results,F,options)

  
    [box,model_exist] = view_model(model,options);

    util = MSAutils;
    [opt, ~, coordinate_system] = util.get_direction(model);
    Qstr = model.output_config.view.on_model{2};
    
    strbox{1} = 'VISUALIZATION:';
    strbox{2} = ' ';
    strbox{3} = ['ANALYSIS: ' model.analysis_type];
    c = 4;
    if options.deformed
        strbox{c} = ['DATA ON MODEL: ' opt ' ' Qstr];
        strbox{c+1} = ['COORDINATE: ' coordinate_system];
        strbox{c+2} = 'SHOW DEFORMED: true';
        c = c + 3;

        %options = struct('color_nodes',[0.5 0.5 0.5],'color_elements',[0,0,0],'color_support',[0.5 0.5 0.5],'holdon'...
        %            ,true,'node_num',out.view.id,'elem_num',out.view.id,'alpha_col',0.3);

        %[box,model_exist] = view_model(model,options);

        options_def = struct('color_nodes',[0 0 0],'color_elements',[0,0,0],'color_support','k','holdon'...
                    ,true,'node_num',false,'elem_num',false,'alpha_col',0.3);

        [~,~] = view_model_deformed(model,results,options_def);

    else
        strbox{c} = 'SHOW DEFORMED: false';
        c = c + 1;

    end


 
    if model_exist

        if isempty(F)
            if options.view_forces 
                color = 'b';
                view_forces(results,results.F_ext,box,options.numerical,color)
            end
    
            if options.view_reactions
                color = 'r';
                view_forces(results,results.R,box,options.numerical,color)
            end

            
            if options.diagrams_on_model

                kind = options.diagram_kind;
                strbox{c} = ['DIAGRAM: ' kind];
                c = c + 1;
                diagrams(results,kind);

            end

        else
            if numel(F)~=numel(results.F_ext)
                error('MSAtool: Cannot plot the model with the given F vector.')
            else
                color = 'b';
                view_forces(results,F,box,options.numerical,color)
            end
        end
        
        if options.box
            hold on 
            dim = [0.01 0.5 .2 .29];
            annotation('textbox',dim,'String',strbox)
        end

        if options.print
            print(options.fig_name,'-dpng')
        end

    else
        disp('MSAtool: no model is shown.')
    end

end


function msa_viewer(model,results)

    util = MSAutils;
    out = model.output_config;
    
    [opt, ~, coordinate_system] = util.get_direction(model);
    Qstr = model.output_config.view.on_model{2};

    strbox{1} = 'VISUALIZATION:';
    strbox{2} = ' ';
    strbox{3} = ['ANALYSIS: ' model.analysis_type];
    c = 4;
    if out.visualization
        
        if out.view.deformed
            strbox{c} = ['DATA ON MODEL: ' opt ' ' Qstr];
   
            strbox{c+1} = ['COORDINATE: ' coordinate_system];
            strbox{c+2} = 'SHOW DEFORMED: yes';
            c = c + 3;

            options = struct('color_nodes',[0.5 0.5 0.5],'color_elements',[0,0,0],'color_support',[0.5 0.5 0.5],'holdon'...
                        ,true,'node_num',out.view.id,'elem_num',out.view.id,'alpha_col',0.3);

            [box,model_exist] = view_model(model,options);

            options = struct('color_nodes',[0 0 0],'color_elements',[0,0,0],'color_support','k','holdon'...
                        ,true,'node_num',out.view.id,'elem_num',out.view.id,'alpha_col',0.3);

            [~,~] = view_model_deformed(model,results,options);
            %[box,model_exist] = view_model_deformed(results,options);
            %error('MSAtool: deformed visualization not implemented.');
            %model_exist = false;
        else
            strbox{c} = 'SHOW DEFORMED: false';
            c = c + 1;

            options = struct('color_nodes',[0 0 0],'color_elements',[0,0,0],'color_support','k','holdon'...
                        ,true,'node_num',out.view.id,'elem_num',out.view.id,'alpha_col',1);

            [box,model_exist] = view_model(model,options);
        end
        
        if model_exist
            if out.view.forces 
                color = 'b';
                view_forces(results,results.F_ext,box,out.view.numerical,color)
            end
    
            if out.view.reactions
                color = 'r';
                view_forces(results,results.R,box,out.view.numerical,color)
            end

            if out.diagrams.on_model{1}
                kind = out.diagrams.on_model{2};

                strbox{c} = ['DIAGRAM: ' kind];
                c = c + 1;

                diagrams(results,kind);
            end

            hold on 
            dim = [0.01 0.5 .2 .29];
            annotation('textbox',dim,'String',strbox)
        else
            disp('MSAtool: no model to show.')
        end

    end
end

function [box,flag] = view_model_deformed(model,results,options)
    
    color_nodes = options.color_nodes;
    color_elements = options.color_elements;
    color_support = options.color_support;
%     holdon = options.holdon;
%     node_num = options.node_num;
%     elem_num = options.elem_num;
%     alpha_col = options.alpha_col;
%     out = model.output_config;

    %util = MSAutils;

    nelem = results.number_elements;
    box = [];

    L = zeros(nelem,1);
    Ulin = zeros(nelem,4);
    %Uang = zeros(nelem,2);
    for i=1:nelem
        L(i) = results.element(i).L;
        %Ulin(i,:) =  results.element(i).U([1 2 4 5])';
        Ulin(i,:) =  results.element(i).U([1 2 3 4])';
        %Uang(i,:) =  results.element(i).U([3 6])';
    end

    l_max = max(L);
    %u_max = max(abs(Ulin(:)));
    %r_max = max(abs(Uang(:)));

    %rnode = l_max/80;
    %rhinge = l_max/60;
    %rh = rnode+rhinge;

    %scale = 0.1*l_max/u_max;

    %[opt, opt_type, coordinate_system] = util.get_direction(model.output_config.view.on_model{2});

    rnode = l_max/80;
    rhinge = l_max/60;
    rh = rnode+rhinge;


    node_ = results.node;
    for i=1:results.number_elements

        xg = results.element(i).x_deformed;
        yg = results.element(i).y_deformed;
     
        hold on
        plot(xg,yg,'Color',color_elements,'Linewidth',2)

%         if results.element(i).node_1.hinge==1
%             plot_node(xg(1),yg(1),rhinge,[1 1 1]);
%         else
%             plot_node(xg(end),yg(end),rnode,color_nodes);
%         end
% 
%         if results.element(i).node_2.hinge==1
%             plot_node(xg(1),yg(1),rhinge,[1 1 1]);
%         else
%             plot_node(xg(end),yg(end),rnode,color_nodes);
%         end
      

        if model.element(i).hinge(1)==1
            if model.element(i).node_1.hinge==0

                %t = model.element(i).angle;
                t = atan2d((yg(end)-yg(end-1)),(xg(end)-xg(end-1)));
                xh = rh*cosd(t) + xg(1);
                yh = rh*sind(t) + yg(1);

                plot_node(xg(1),yg(1),rnode,color_nodes);
                plot_node(xh,yh,rhinge,[1 1 1]);

            else
                plot_node(xg(1),yg(1),rhinge,[1 1 1]);
            end

        else

            if model.element(i).node_1.hinge==0
                plot_node(xg(1),yg(1),rnode,color_nodes);
            else
                plot_node(xg(1),yg(1),rhinge,[1 1 1]);
            end

        end

        if model.element(i).hinge(2)==1
            if model.element(i).node_2.hinge==0

                %model.element(i).angle
                t0 = atan2d((yg(end)-yg(end-1)),(xg(end)-xg(end-1)));
                t = 180+t0;
                xh = rh*cosd(t) + xg(end);
                yh = rh*sind(t) + yg(end);

                plot_node(xg(end),yg(end),rnode,color_nodes);
                plot_node(xh,yh,rhinge,[1 1 1]);

            else
                plot_node(xg(end),yg(end),rhinge,[1 1 1]);
            end

        else

            if model.element(i).node_2.hinge==0
                plot_node(xg(end),yg(end),rnode,color_nodes);
            else
                plot_node(xg(end),yg(end),rhinge,[1 1 1]);
            end

        end

        node_(results.element(i).node_1.id).x = xg(1);
        node_(results.element(i).node_1.id).y = yg(1);
        node_(results.element(i).node_2.id).x = xg(end);
        node_(results.element(i).node_2.id).y = yg(end);

        x(results.element(i).node_1.id) = xg(1);
        x(results.element(i).node_2.id) = xg(end);
        y(results.element(i).node_1.id) = yg(1);
        y(results.element(i).node_2.id) = yg(end);

        n = length(xg);
 
        if model.output_config.view.on_model{1}
            
            %coordinate_system
            c = results.element(i).color_results;
    
            %Flim = results.element(i).color_results_limits;
          
            %Q = results.element(i).Q;

            if results.element(i).rigid
                
                patch('faces',[1:n-1;2:n]','vertices',[xg; yg]',...
                      'faceVertexCData',c,...
                      'edgecolor','flat',...
                      'linewidth',6);
            else
                patch('faces',[1:n-1;2:n]','vertices',[xg; yg]',...
                  'faceVertexCData',c,...
                  'edgecolor','flat',...
                  'linewidth',2);
            end
            
        end

    end

    if model.output_config.view.on_model{1}
        cbar = results.colorbar;
        hcb=colorbar;
        hcb.Title.String = model.output_config.view.on_model{2};
        %colorbar;
        colormap jet;
    
        %set(gca,'CLim',[Fmin Fmax]);
        if cbar(1) ~= cbar(2)    
            set(gca,'CLim',[cbar(1) cbar(2)]);
        else
            
            if cbar(1)==0
                set(gca,'CLim',[-1 1]);
            else
                set(gca,'CLim',[-cbar(1) cbar(2)]);
            end
        end
        %annotation('textbox', [0, 0.5, 0, 0], 'string', coordinate_system)
        %annotation('textbox', [0, 0.5, 0, 0], 'string', opt)
         %dim = [0.01 0.5 .2 .2];
         %str = {[coordinate_system ' ' opt]};
         %str = {'aaaaa'};
        %annotation('textbox',dim,'String',str)
    end

    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    mxx = max([maxx,maxy]);
    mix = max([minx,miny]);

    maxx = mxx;
    minx = mix;
    %maxy = mxx;
    %miny = mix;

    % Determine the size of each dimension.
    Dx = maxx - minx;
    %Dy = maxy - miny;

    for i=1:numel(model.support)
        node1 = node_(model.support(i).node_id);
        plot_support(model.support(i),node1,Dx,color_support);
    end

    flag = true;

end

function [box,flag] = view_model(model,options)

    dim = 2;
    color_nodes = options.color_nodes;
    color_elements = options.color_elements;
    color_support = options.color_support;
    holdon = options.holdon;
    node_num = options.node_num;
    elem_num = options.elem_num;
    alpha_col = options.alpha_col;

    x = zeros(1,model.number_nodes);
    y = zeros(1,model.number_nodes);
    for i=1:model.number_nodes
        x(i) = model.node(i).x;
        y(i) = model.node(i).y;
    end
    
    % Find the minimum and maximum coordinates.
    minx = min(x);
    maxx = max(x);
    miny = min(y);
    maxy = max(y);
    
    mxx = max([maxx,maxy]);
    mix = max([minx,miny]);

    maxx = mxx;
    minx = mix;
    maxy = mxx;
    miny = mix;

    % Determine the size of each dimension.
    Dx = maxx - minx;
    Dy = maxy - miny;

    box = [Dx Dy 0];
    
    util = MSAutils;
    str2rgb = util.str2rgb;

    fid = {};
    for i=1:model.number_nodes
        x(i) = model.node(i).x;
        y(i) = model.node(i).y;
    end
   
    if ~holdon
        %disp('View model: OFF')
        %h=figure('color',[1,1,1],'visible','off');
        %saveas(h,'fig.png');
        figure('color',[1,1,1],'visible','off');
        hold on
    else
        hold on
    end
    
    for i=1:model.number_elements
        Li(i) = model.element(i).L;
    end

    Lmax = max(Li);
    rnode = Lmax/80;
    rhinge = Lmax/60;
    rh = rnode+rhinge;

    for i=1:model.number_elements
        
        L = model.element(i).L;
        x1=model.element(i).node_1.x;
        y1=model.element(i).node_1.y;
        x2=model.element(i).node_2.x;
        y2=model.element(i).node_2.y;

        %rgbVector = str2rgb(color_elements);
        rgb_alpha = [color_elements,alpha_col];
        if model.element(i).rigid
            rgb_alpha = [[1 0.5 0.2],max(alpha_col,0.7)];
            plot([x1 x2],[y1 y2],'Color',rgb_alpha,'linewidth',6);
        else
            plot([x1 x2],[y1 y2],'Color',rgb_alpha,'linewidth',2);
        end

        if elem_num
            str = num2str(model.element(i).id);
            text((x1+x2)/2, (y1+y2)/2+Dy/20, ['$\raisebox{.5pt}{\fbox{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex', 'FontSize',14);
        end


        if model.element(i).hinge(1)==1
            if model.element(i).node_1.hinge==0
                x = model.element(i).node_1.x;
                y = model.element(i).node_1.y;

                t = model.element(i).angle;
                xh = rh*cosd(t) + x;
                yh = rh*sind(t) + y;

                plot_node(xh,yh,rhinge,[1 1 1]);
            end
        end

        if model.element(i).hinge(2)==1
            if model.element(i).node_2.hinge == 0
                x = model.element(i).node_2.x;
                y = model.element(i).node_2.y;

                t = 180+model.element(i).angle;
                xh = rh*cosd(t) + x;
                yh = rh*sind(t) + y;

                plot_node(xh,yh,rhinge,[1 1 1]);
            end
        end
       
        
    
    end

    for i=1:model.number_nodes
    
        %plot_node(model.node(i),dim,color_nodes);
        x = model.node(i).x;
        y = model.node(i).y;

        if model.node(i).hinge==1
            plot_node(x,y,rhinge,[1 1 1]);
        else
            plot_node(x,y,rnode,color_nodes);
        end

        if node_num
            str = num2str(model.node(i).id);
            x = model.node(i).x;
            y = model.node(i).y;

            text(x, y+Dy/20, ['$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {' str '}}}$'], 'Interpreter', 'latex','Color','m', 'FontSize',16);

        end
    
    end
 
    
    for i=1:numel(model.support)
         node = model.node(model.support(i).node_id);
         plot_support(model.support(i),node,Dx,color_support);
    end
     
    xlim([minx - Dx/3 maxx + Dx/3]);
    if maxy ~=miny
        ylim([miny - Dy/3 maxy + Dy/3]);
    else
        ylim([miny - Dx/3 maxy + Dx/3]);
    end
     
    xlabel('X');
    ylabel('Y');

    axis equal;
    axis off;

    flag = true;
    %if printfig
    %    print(fig_model,'-dpng')
    %end


end

function view_forces(results,F,box,numerical,color)

    for i=1:results.number_nodes
            
            %force = F(results.node(i).dof);
            %force(dofn) = F(dofg);
            position = [results.node(i).x results.node(i).y 0];
            direc = results.node(i).directions;
            force = zeros(3,1);
            force(direc) = F(results.node(i).dof);
 
            plot_forces(force,position,box,2,color,numerical);
           
    end
    
   
    axis equal;
    axis off;

    %if printfig
    %    print(fig_model,'-dpng')
    %end


end


function diagrams(results,kind)

    util = MSAutils;

    if strcmp(results.analysis_type,'frame')
        if strcmp(kind,'axial')
            diagram_option = 1;
        elseif strcmp(kind,'shear')
            diagram_option = 2;
        elseif strcmp(kind,'bending')
            diagram_option = 3;
        else
            error('MSAtool: not recognized type of internal force diagram. The existing options are: axial, shear, bending.');
        end

    elseif strcmp(results.analysis_type,'beam')

        if strcmp(kind,'shear')
            diagram_option = 1;
        elseif strcmp(kind,'bending')
            diagram_option = 2;
        else
            error('MSAtool: not recognized type of internal force diagram. The existing options are: axial, shear, bending.');
        end

    elseif strcmp(results.analysis_type,'truss') || strcmp(results.analysis_type,'bar')

        if strcmp(kind,'axial')
            diagram_option = 1;
        else
            error('MSAtool: not recognized type of internal force diagram. The existing options are: axial, shear, bending.');
        end

    end


    for i=1:results.number_elements
     
        %direc = {'x';'y';'xy'};
        P = results.element(i).P;
        %P1 = P(1:3);
        %P2 = P(4:6);
        nforce = numel(results.element(i).internal_force(:,1));
        c = 1;
        for j=1:nforce
            iforce(i,j,:) = results.element(i).internal_force(j,:);
            c = c + 1;

        end

        Lall(i) = results.element(i).L;

    end

    fo = iforce(:,diagram_option,:);
    maxf = max(abs(fo(:)));
    maxL = max(Lall);
    scale = 0.1*maxL/maxf;
    dx = maxL/10; %30 before

    for i=1:results.number_elements
        
        x1 = results.element(i).node_1.x;
        x2 = results.element(i).node_2.x;
        y1 = results.element(i).node_1.y;
        y2 = results.element(i).node_2.y;

        t = results.element(i).angle;
        R = [cosd(t) -sind(t);sind(t) cosd(t)];

        Y0 = iforce(i,diagram_option,:);
        Y = Y0*scale;
        X1 = [0;Y(1)];
        X2 = [results.element(i).L;Y(2)];
        
        X1p = R*X1 + [x1;y1];
        X2p = R*X2 + [x1;y1];

        sgn = sign([Y0(1) Y0(2)]);
          
        if sgn(1) == sgn(2)

            hold on
            %plot([X1p(1) X2p(1)],[X1p(2) X2p(2)])
            %plot([x1 x2],[y1 y2]);
            xx = [x1 x2 X2p(1) X1p(1)];
            yy = [y1 y2 X2p(2) X1p(2)];
            pgon = polyshape(xx,yy);

            if sgn(1) < 0
                plot(pgon,'EdgeColor',[0 0 0],'FaceColor',[1 0.5 0],'FaceAlpha',0.1);
            else
                plot(pgon,'EdgeColor',[0 0 1],'FaceColor',[0 0.5 1],'FaceAlpha',0.1);
            end

        else

            p1i = [x1 y1];
            p1f = [x2 y2];
            p2i = [X1p(1) X1p(2)];
            p2f = [X2p(1) X2p(2)];
            a1 = (p1f(2) - p1i(2))/(p1f(1) - p1i(1));
            a2 = (p2f(2) - p2i(2))/(p2f(1) - p2i(1));
            b1 = -a1*p1i(1) + p1i(2);
            b2 = -a2*p2i(1) + p2i(2);

    
            if isinf(a1) && ~isinf(a2)
    
                xc = p1i(1);
                yc = a2*xc + b2;
                
            elseif ~isinf(a1) && isinf(a2)
    
                xc = p2i(1);
                yc = a1*xc + b1;
    
            elseif isinf(a1) && isinf(a2)
    
                xc = inf;
                yc = inf;
    
            else
    
                xc = -(b2 - b1)/(a2 - a1+0.0);
                yc = a1*xc + b1;
    
            end

            hold on
            %plot([X1p(1) xc],[X1p(2) yc])
            %plot([x1 xc],[y1 yc]);
            xx = [x1 xc xc X1p(1)];
            yy = [y1 yc yc X1p(2)];
  
            pgon = polyshape(xx,yy);

            if sgn(1) < 0
                plot(pgon,'FaceColor',[1 0.5 0],'FaceAlpha',0.1);
            else
                plot(pgon,'FaceColor',[0 0.5 1],'FaceAlpha',0.1);
            end
            % 
            hold on
            %plot([xc X2p(1)],[yc X2p(2)])
            %plot([xc x2],[yc y2]);
            xx = [xc x2 X2p(1) xc];
            yy = [yc y2 X2p(2) yc];
            pgon = polyshape(xx,yy);

            if sgn(2) < 0
                plot(pgon,'FaceColor',[1 0.5 0],'FaceAlpha',0.1);
            else
                plot(pgon,'FaceColor',[0 0.5 1],'FaceAlpha',0.1);
            end

        end


        if X1p(1) < X2p(1)
            mx = [-1 1]*dx;
        else
            mx = [1 -1]*dx;
        end

        if X1p(2) > 0
            my1 = dx;
        else
            my1 = -dx;
        end

        if X2p(2) > 0
            my2 = dx;
        else
            my2 = -dx;
        end
        
        text(X1p(1)+mx(1),X1p(2)+my1,num2str(round(Y0(1),2)),'Color','r','FontSize',12);
        text(X2p(1)+mx(2),X2p(2)+my2,num2str(round(Y0(2),2)),'Color','r','FontSize',12);
  
    end

%     if diagram_option==1
%         tx = 'Axial force';
%     elseif diagram_option==2
%         tx = 'Shear force';
%     elseif diagram_option==3
%         tx = 'Bending moment';
%     end
% 
%     title(tx,'FontSize',18)

end


function plot_node(x,y,r,color_nodes)

    %sc=scatter(x,y,'MarkerEdgeColor',[0 0 0],...
    %          'MarkerFaceColor',color_nodes, 'LineWidth',1.5);

    %sc.SizeData = nodesize;
    %viscircles([x,y],1,'Color',[0 0 0]);

    %// radius
    %r = 1;
    
    %// center
    c = [x y];
    
    pos = [c-r 2*r 2*r];
    %r = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor', 'red',
    %'Edgecolor','none','FaceAlpha',1);% only in 2024
    
    r = rectangle('Position',pos,'Curvature',[1 1], 'FaceColor',color_nodes, 'Edgecolor',[0 0 0]);

    hold on

end

function plot_forces(force,position,Dxyz,dim,color,numerical)

    util = MSAutils;
    circular_arrow = util.circular_arrow;
  
    Dx = Dxyz(1);
    Dy = Dxyz(2);
    Dz = Dxyz(3);

    x = position(1);
    y = position(2);
    z = position(3);

    force(isnan(force))=0;

    %idx = [1 2 6];
    F(1).force = [force(1) 0];
    F(2).force = [0 force(2)];
    F(3).force = force(3);
    F(1).type = 'F';
    F(2).type = 'F';
    F(3).type = 'M';
    num_force = 3;

    for i=1:num_force
        
        f = F(i).force;
        fnorm = norm(f);
        if fnorm==0
            fnorm=1;
        end
        fn = f/fnorm;
        
        if strcmp(F(i).type,'F')

            if dim==1

                if sum(fn) ~= 0
                    p1 = [x-(Dx/5)*fn 0];
                    p2 = [x-Dx/50*fn  0];
                    
                    dp = p2-p1;    
                    quiver(p1(1),p1(2),dp(1),dp(2),0,color,'LineWidth',4,'MaxHeadSize',1);  
                end

                if sum(fn) ~= 0 && numerical
                    text(p1(1)+fn(1)*Dx/50,p1(2)+Dx/50, sprintf(num2str(sum(abs(f))),p1),'Color',color,'FontSize',14);
                end

            elseif dim==2

                if sum(fn) ~= 0
                    p2 = [x-(Dx/50)*fn(1)  y-(Dy/50)*fn(2)];
                    p1 = [x-(Dx/5)*fn(1)  y-(Dy/5)*fn(2)];
        
                    dp = p2-p1;    
                    quiver(p1(1),p1(2),dp(1),dp(2),0,color,'LineWidth',1,'MaxHeadSize',1); 
               
                end

                if sum(fn) ~= 0 && numerical
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

            end 

            if f ~= 0 && numerical
                
                p1 = [x+Dx/12 y+Dy/12];
                text(p1(1),p1(2), sprintf(num2str(sum(abs(f))),p1),'Color',color,'FontSize',14)
            end 
        else

        end

    end

end

% function plot_support_12D(support,color_support,Dx,dim)
% 
%     hold on
%     draw_fixed_support_12D(Dx,support,color_support);
% 
% end

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

function plot_support(support,node,Dx,color_support)
    

    Fcol = color_support;
   
    constr = support.constraint';
    theta = support.angle;
    
    R = [cosd(theta) -sind(theta);sind(theta) cosd(theta)];

    side_len = Dx/15;
    h = side_len*sin(pi/3)/1.3;
    h0 = side_len*sin(pi/3);
    centerpoint = [node.x ,node.y-h];
    p = nsidedpoly(3, 'Center', centerpoint, 'SideLength', side_len);

    refpoint = [node.x ,node.y];
    p = rotate(p,theta,refpoint);


    if isequal(constr,[0 1 1])

        hold on

        px = [node.x-(h0/1.5), node.x+(h0/1.5)];
        py = [node.y-0.2*h node.y-0.2*h];
        pyy = [node.y-0.7*h node.y-0.7*h];

        c = [node.x;node.y];
        p1 = [px(1);py(1)];
        p2 = [px(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr1 = [p1r(1) p2r(1)];
        pyr1 = [p1r(2) p2r(2)];

        p1 = [px(1);pyy(1)];
        p2 = [px(2);pyy(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr2 = [p1r(1) p2r(1)];
        pyr2 = [p1r(2) p2r(2)];

        plot(pxr1,pyr1,'Color',color_support,'LineWidth',2);
        plot(pxr2,pyr2,'Color',color_support,'LineWidth',2);
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 0 0])

        hold on

        px = [node.x-0.2*h node.x-0.2*h];
        pxx = [node.x+0.2*h node.x+0.2*h];
        py = [node.y-(h0/1.5) node.y+(h0/1.5)];


        c = [node.x;node.y];
        p1 = [px(1);py(1)];
        p2 = [px(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr1 = [p1r(1) p2r(1)];
        pyr1 = [p1r(2) p2r(2)];

        p1 = [pxx(1);py(1)];
        p2 = [pxx(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr2 = [p1r(1) p2r(1)];
        pyr2 = [p1r(2) p2r(2)];

        plot(pxr1,pyr1,'Color',color_support,'LineWidth',2);
        plot(pxr2,pyr2,'Color',color_support,'LineWidth',2);
        pg.FaceColor = Fcol;

     elseif isequal(constr,[1 0 1])

        hold on

        px = [node.x-0.2*h node.x-0.2*h];
        pxx = [node.x+0.2*h node.x+0.2*h];
        pxc = [node.x node.x];
        py = [node.y-(h0/1.5) node.y+(h0/1.5)];

        %px = [node.x-(h0/1.5), node.x+(h0/1.5)];
        %py = [node.y-0.2*h node.y-0.2*h];
        %pyy = [node.y-0.7*h node.y-0.7*h];

        c = [node.x;node.y];
        p1 = [px(1);py(1)];
        p2 = [px(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr1 = [p1r(1) p2r(1)];
        pyr1 = [p1r(2) p2r(2)];

        p1 = [pxx(1);py(1)];
        p2 = [pxx(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr2 = [p1r(1) p2r(1)];
        pyr2 = [p1r(2) p2r(2)];

        p1 = [pxc(1);py(1)];
        p2 = [pxc(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxrc = [p1r(1) p2r(1)];
        pyrc = [p1r(2) p2r(2)];

        plot(pxr1,pyr1,'Color',color_support,'LineWidth',2);
        plot(pxr2,pyr2,'Color',color_support,'LineWidth',2);
        plot(pxrc,pyrc,'Color',color_support,'LineWidth',2);
        pg.FaceColor = Fcol;

    elseif isequal(constr,[0 1 0])


        hold on

        px = [node.x-(h0/1.5), node.x+(h0/1.5)];
        py = [node.y-1.8*h node.y-1.8*h];

        c = [node.x;node.y];
        p1 = [px(1);py(1)];
        p2 = [px(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr1 = [p1r(1) p2r(1)];
        pyr1 = [p1r(2) p2r(2)];

        
    
        pg = plot(p);
        %plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-1.8*h node.y-1.8*h],'Color',color_support,'LineWidth',2)
        plot(pxr1,pyr1,'Color',color_support,'LineWidth',2);
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 1 0])

        pg = plot(p);
        pg.FaceColor = Fcol;

    elseif isequal(constr,[1 1 1])

        hold on
        
        
        px = [node.x-(h0/1.5), node.x+(h0/1.5)];
        py = [node.y-0.2*h node.y-0.2*h];

        c = [node.x;node.y];
        p1 = [px(1);py(1)];
        p2 = [px(2);py(2)];
        p1r = R*(p1-c) + c;
        p2r = R*(p2-c) + c;
        pxr1 = [p1r(1) p2r(1)];
        pyr1 = [p1r(2) p2r(2)];

        plot(pxr1,pyr1,'Color',color_support,'LineWidth',2);
        %plot([node.x-(h0/1.5), node.x+(h0/1.5)],[node.y-0.5*h node.y-0.5*h],'Color',color_support,'LineWidth',2)
        %pg.FaceColor = Fcol;

    else

        error('MSAtool error: inconsistent supports!')

    end
    
end

