function opt_vars = extract_opt_vars(x_opt, params, theta_0, d_end_goal)

    N_pts = params.N_pts;
    N_links = params.N_links;
    L_0 = params.L_0;

    idx_x = 0; % index into x_opt
    opt_vars = struct();

    % loop over links
    for link = 1:N_links
        S_Px_opt  = zeros(1, N_pts);
        S_Py_opt  = zeros(1, N_pts);
        S_Dx_opt  = zeros(1, N_pts);
        S_Dy_opt  = zeros(1, N_pts);
        alpha_P_opt   = zeros(1, N_pts);
        alpha_D_opt   = zeros(1, N_pts);
        theta_opt  = zeros(1, N_pts);

        for k = 1:N_pts
            % extract decision variables in correct order
            idx_x = idx_x + 1;
            S_Px_opt(k) = x_opt(idx_x);
            idx_x = idx_x + 1;
            S_Py_opt(k) = x_opt(idx_x);
            idx_x = idx_x + 1;
            S_Dx_opt(k) = x_opt(idx_x);
            idx_x = idx_x + 1;
            S_Dy_opt(k) = x_opt(idx_x);
            idx_x = idx_x + 1;
            alpha_P_opt(k)  = x_opt(idx_x);
            idx_x = idx_x + 1;
            alpha_D_opt(k)  = x_opt(idx_x);
            idx_x = idx_x + 1;
            theta_opt(k) = x_opt(idx_x);
        end

        % store extracted variables for this link
        opt_vars.S(link).S_Px = S_Px_opt;
        opt_vars.S(link).S_Py = S_Py_opt;
        opt_vars.S(link).S_Dx = S_Dx_opt;
        opt_vars.S(link).S_Dy = S_Dy_opt;
        opt_vars.S(link).alpha_P  = alpha_P_opt;
        opt_vars.S(link).alpha_D  = alpha_D_opt;
        opt_vars.S(link).theta = theta_opt;
    end

    % reconstruct d_end and pulley points based on decision variables
    dendx_opt   = zeros(1, N_pts);
    dendy_opt   = zeros(1, N_pts);
    dend_theta  = zeros(1, N_pts);

    for link = 1:N_links

        L = params.link(link).L;
        S_Px_opt = opt_vars.S(link).S_Px;
        S_Py_opt = opt_vars.S(link).S_Py;
        S_Dx_opt = opt_vars.S(link).S_Dx;
        S_Dy_opt = opt_vars.S(link).S_Dy;
        theta_opt = opt_vars.S(link).theta;

        for k = 1:N_pts
            Rtheta = rotate_2D(theta_opt(k));
            S_DRef = [S_Px_opt(k); S_Py_opt(k)] - Rtheta * [S_Dx_opt(k); S_Dy_opt(k)];
            old_angle = dend_theta(k);
            new_angle = old_angle + theta_opt(k);
            temp = rotate_2D(old_angle) * S_DRef;
            dendx_opt(k) = dendx_opt(k) + temp(1) + L * cos(new_angle);
            dendy_opt(k) = dendy_opt(k) + temp(2) + L * sin(new_angle);
            dend_theta(k) = new_angle;
        end
    end

    opt_vars.dendx = dendx_opt;
    opt_vars.dendy = dendy_opt;

    for link = 1:N_links
        alpha_P_opt  = opt_vars.S(link).alpha_P;
        alpha_D_opt  = opt_vars.S(link).alpha_D;
        theta_opt = opt_vars.S(link).theta;
        S_Px_opt = opt_vars.S(link).S_Px;
        S_Py_opt = opt_vars.S(link).S_Py;
        S_Dx_opt = opt_vars.S(link).S_Dx;
        S_Dy_opt = opt_vars.S(link).S_Dy;
        R  = params.link(link).R;

        if link == 1
            gamma_P = alpha_P_opt - theta_0; 
        else
            gamma_P = opt_vars.S(link-1).theta + alpha_P_opt;
        end
        gamma_D = alpha_D_opt - theta_opt;

        % frame_D written in frame_P
        S_DRef_opt = zeros(2, N_pts);
        for k = 1:N_pts
            Rtheta = rotate_2D(theta_opt(k));
            S_DRef_opt(:, k) = [S_Px_opt(k); S_Py_opt(k)] - Rtheta * [S_Dx_opt(k); S_Dy_opt(k)];
        end
        
        % calculate pulley line constants
        c_D_opt = zeros(1, N_pts);
        for k = 1:N_pts
            c_D_opt(k) = [-sin(alpha_D_opt(k)), cos(alpha_D_opt(k))] * S_DRef_opt(:, k) + params.link(link).phi_D*R;
        end

        c_P_opt = zeros(1, N_pts);
        if link == 1
            for k = 1:N_pts
                c_P_opt(k) = L_0 * sin(alpha_P_opt(k)) + params.link(1).phi_D*R;
            end
        else
            S_DRef_prev = opt_vars.S(link-1).S_DRef;
            for k = 1:N_pts
                c_P_opt(k) = [sin(gamma_P(k)), -cos(gamma_P(k))] * S_DRef_prev(:, k) + params.link(link).phi_P*R;
            end
        end

        % reconstructed pulley shapes in local frames
        P_P = zeros(2, N_pts-1);
        P_D = zeros(2, N_pts-1);
        for k = 2:N_pts
            P_P(:,k-1) = [sin(gamma_P(k-1)), -cos(gamma_P(k-1)); sin(gamma_P(k)), -cos(gamma_P(k))] \ [c_P_opt(k-1); c_P_opt(k)];
            P_D(:,k-1) = [sin(gamma_D(k-1)), -cos(gamma_D(k-1)); sin(gamma_D(k)), -cos(gamma_D(k))] \ [c_D_opt(k-1); c_D_opt(k)];
        end

        opt_vars.S(link).S_DRef = S_DRef_opt;
        opt_vars.S(link).P_P = P_P;
        opt_vars.S(link).P_D = P_D;

    end

    opt_vars.d_end_goal = d_end_goal;
    opt_vars.theta_0 = theta_0;
end
