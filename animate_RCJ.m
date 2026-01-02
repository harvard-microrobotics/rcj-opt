function animate_RCJ(opt_vars, params, filename)

    N_links = params.N_links;
    N_pts   = params.N_pts;
    Colors = {[205 38 39]/255, [55 76 158]/255, [250 168 54]/255, [205 38 39]/255, [55 76 158]/255};
    
    % figure setup
    fig = figure;
    hold on 
    axis equal
    box on
    set(gca, 'LineWidth', 2);
    xlabel('x')
    ylabel('y')
    xmin = -6;
    xmax = 8;
    ymin = -6;  
    ymax = 6;
    axis([xmin xmax ymin ymax]);
    
    % optimized end effector path
    plot(opt_vars.dendx, opt_vars.dendy, '--', 'Color','black','LineWidth',2,'DisplayName','End Effector Path');

    % store animation handles
    S_P_curve = gobjects(N_links,1); % surfaces
    S_D_curve = gobjects(N_links,1);
    P_P_curve = gobjects(N_links,1); % pulleys
    P_D_curve = gobjects(N_links,1);
    bearing_P_curve = gobjects(N_links,1); % idler bearings
    bearing_D_curve = gobjects(N_links,1);

    % initialize plots
    for link = 1:N_links
        S_P_curve(link) = plot(nan, nan, 'Color', Colors{link}, 'LineWidth', 2);
        S_D_curve(link) = plot(nan, nan, 'Color', Colors{link+1}, 'LineWidth', 2);
        P_P_curve(link) = plot(nan, nan, 'Color', Colors{link}, 'LineWidth', 2);
        P_D_curve(link) = plot(nan, nan, 'Color', Colors{link+1}, 'LineWidth', 2);
        bearing_D_curve(link) = plot(nan, nan, 'k', 'LineWidth', 1);
        bearing_P_curve(link) = plot(nan, nan, 'k', 'LineWidth', 1);
    end
    end_point = scatter(nan, nan, 50, 'black', 'filled');

    % store video frames
    frames(N_pts) = struct('cdata',[],'colormap',[]);
    frame_D_prev = [0;0;0];
    for i = 1:N_pts % frame
        
        % frames are the same for link 1 and start at the origin
        frame_P = [0; 0; 0];
        frame_D = [0; 0; 0];

        for link = 1:N_links
            
            % shorthand
            S_Px = opt_vars.S(link).S_Px;
            S_Py = opt_vars.S(link).S_Py;
            S_Dx = opt_vars.S(link).S_Dx;
            S_Dy = opt_vars.S(link).S_Dy;
            P_P = opt_vars.S(link).P_P;
            P_D = opt_vars.S(link).P_D;
            theta_0 = opt_vars.theta_0;
            theta = opt_vars.S(link).theta;
            
            L_0 = params.L_0;
            L = params.link(link).L;
            R = params.link(link).R;

            % bearing circle points (untransformed)
            s  = linspace(0,2*pi,100);
            bx = R * cos(s);
            by = R * sin(s);
            bearing_pts = [bx; by];

            % frame rotate_2Dations and trnslations
            Rk_P2 = rotate_2D(theta(i));
            T_P2  = [S_Px(i); S_Py(i)] - rotate_2D(theta(i)) * [S_Dx(i); S_Dy(i)];

            if link == 1
                Rk_P1 = rotate_2D(theta_0(i));
                T_P1  = [-L_0; 0];
            else
                Rk_P1 = rotate_2D(opt_vars.S(link-1).theta(i));
                T_P1  = frame_D(1:2);
            end

            % transform surfaces and pulleys to global frame
            S_P_global = zeros(2, N_pts);
            S_D_global = zeros(2, N_pts);
            P_P_global = zeros(2, N_pts-1);
            P_D_global = zeros(2, N_pts-1);

            for k = 1:N_pts
                S_P_global(:,k) = rigid_tf_2D([S_Px(k); S_Py(k)], frame_D);
                S_D_global(:,k) = rigid_tf_2D(Rk_P2 * [S_Dx(k); S_Dy(k)] + T_P2, frame_D);
            end

            for k = 1:N_pts-1
                if link == 1
                    P_P_global(:,k) = Rk_P1 * P_P(:,k) + T_P1;
                else
                    P_P_global(:,k) = rigid_tf_2D(P_P(:,k), frame_D_prev);
                end    
                P_D_global(:,k) = rigid_tf_2D(Rk_P2 * P_D(:,k) + T_P2, frame_D);
            end

            set(P_P_curve(link), 'XData', P_P_global(1,:), 'YData', P_P_global(2,:));
            set(P_D_curve(link), 'XData', P_D_global(1,:), 'YData', P_D_global(2,:));
            set(S_P_curve(link), 'XData', S_P_global(1,:), 'YData', S_P_global(2,:));
            set(S_D_curve(link), 'XData', S_D_global(1,:), 'YData', S_D_global(2,:));

            % transform bearing 1 to global frame
            bearing_global = rigid_tf_2D(bearing_pts, frame_D);
            set(bearing_D_curve(link), 'XData', bearing_global(1,:), 'YData', bearing_global(2,:));

            % advance frames
            frame_D_prev = frame_D;
            frame_P(1:2) = frame_D(1:2) + rotate_2D(frame_D(3)) * T_P2;
            frame_P(3) = frame_P(3) + theta(i);
            frame_D(1:2) = frame_P(1:2) + rotate_2D(frame_P(3)) * [L; 0];
            frame_D(3) = frame_D(3) + theta(i);

            % bearing 2 plot
            bearing_global2 = rigid_tf_2D(bearing_pts, frame_P);
            if link < N_links
                set(bearing_P_curve(link), 'XData', bearing_global2(1,:), 'YData', bearing_global2(2,:));
            end

            if link == N_links
                set(end_point, 'XData', frame_D(1), 'YData', frame_D(2));
            end
        end
        
        % render and store frames
        drawnow;
        frames(i) = getframe(fig);
    end
   
    % save video
    video_title = strcat(filename, '.mp4');
    obj = VideoWriter(video_title, 'MPEG-4');
    obj.Quality = 100;
    obj.FrameRate = 10; 
    open(obj);
    iter = 1:N_pts;
    for i = iter % reverse video
        writeVideo(obj, frames(i));
    end
    close(obj);

end

function X_transformed = rigid_tf_2D(X, transform_vec)
    % rotate_2Date and translate X according to transform_vec
    X_transformed = rotate_2D(transform_vec(3)) * X + transform_vec(1:2);
end