function params = generate_params(N_links, N_pts)
    % creates a default params structure for the RCJ optimization with constraint bounds and parameters
    
    params.L_0 = 3.0;  % distance from input pulley to first idler bearing
    params.N_links = N_links;
    params.N_pts = N_pts;
    params.eta_0 = -1; % 1 for CCW -1 for CW rotate_2Dation of input pulley
    params.omega = pi; % total rotate_2Dation of input pulley (rad)

    % constraint activation flags
    params.use_constraints = struct();
    params.use_constraints.S_P_size = true;
    params.use_constraints.S_D_size = true;
    params.use_constraints.S_tangency = true;
    params.use_constraints.S_arc_length = true;
    params.use_constraints.S_step_size = true;
    params.use_constraints.S_P_convexity = true;
    params.use_constraints.S_D_convexity = true;
    params.use_constraints.r_P_length = true;
    params.use_constraints.r_D_length = true;
    params.use_constraints.P_moment_arm = true; 
    params.use_constraints.P_P_rotate_2Dation = true;
    params.use_constraints.P_D_rotate_2Dation = true;
    params.use_constraints.P_P_size = true;
    params.use_constraints.P_D_size = true;
    params.use_constraints.P_D_intercept = true;
    params.use_constraints.P_P_intercept = true;
    params.use_constraints.L_wire = true;


    for link = 1:N_links
        params.link(link).R = 0.4; % idler bearing radii
        params.link(link).L = 1; % distance between idler bearings, or from distal frame to end effector
        params.link(link).eta = -1; % joint rotate_2Dation for each link 1 for CCW -1 for CW

        if link == 1
            % bearing wrap should be opposite joint direction for favorable
            % force transmission in most cases.
            params.link(link).phi_D = -params.link(link).eta;
            params.link(link).phi_P = params.link(link).phi_D;
        else
            params.link(link).phi_P = params.link(link-1).eta;
            params.link(link).phi_D = -params.link(link).eta;
        end

        % decision variable bounds
        params.link(link).var_bnd.S_Px = [-inf, inf];
        params.link(link).var_bnd.S_Py = [-inf, inf];
        params.link(link).var_bnd.S_Dx = [-inf, inf];
        params.link(link).var_bnd.S_Dy = [-inf, inf];
        params.link(link).var_bnd.alpha_P = [-2*pi, 2*pi];
        params.link(link).var_bnd.alpha_D = [-2*pi, 2*pi];
        params.link(link).var_bnd.theta = [-pi, pi];

        % constraint bounds
        params.link(link).const_bnd.S_P_size = [0.5, 3];
        params.link(link).const_bnd.S_D_size = [0.5, 3];
        params.link(link).const_bnd.S_step_size = [0.001, 0.2];
        params.link(link).const_bnd.P_P_size = [0.05, 0.4];
        params.link(link).const_bnd.P_D_size = [0.05, 0.4];
        params.link(link).const_bnd.r_P_length = [-inf, -0.2];
        params.link(link).const_bnd.r_D_length = [0.2, inf];
        params.link(link).const_bnd.S_tangency = [0, 0];
        params.link(link).const_bnd.S_arc_length = [0, 0];
        params.link(link).const_bnd.P_moment_arm = [0, 0];
        params.link(link).const_bnd.S_P_convexity = [1e-8, inf];
        params.link(link).const_bnd.S_D_convexity = [-inf, -1e-8];

        if link == 1
            params.link(link).const_bnd.P_P_rotate_2Dation = [0, inf];
        else
            params.link(link).const_bnd.P_P_rotate_2Dation = [-inf, 0];
        end

        params.link(link).const_bnd.P_D_rotate_2Dation = [0, inf];
        params.link(link).const_bnd.P_P_intercept = [-0.2, -0.001];
        params.link(link).const_bnd.P_D_intercept = [-0.2, -0.005];
        params.link(link).const_bnd.L_wire = [-0.1/N_pts, 0.1/N_pts];
    end
end

