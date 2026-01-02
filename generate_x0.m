function x0_struct = generate_x0(params)
    % generates structure containing an initial guess for the optimization.
    % initial guess corresponds to circular surfaces and pulleys

    N_pts = params.N_pts;
    N_links = params.N_links;
    omega = params.omega;
    x0_struct = struct();

    for link = 1:N_links
        for k = 1:N_pts

            s_k = (k-1)/(N_pts - 1);
            theta_0 = params.link(link).eta*(-omega/2 + s_k*omega);
            
            x0_struct.S(link).S_Px(k) = cos(theta_0/2);
            x0_struct.S(link).S_Py(k) = sin(theta_0/2);
            x0_struct.S(link).S_Dx(k) = -cos(theta_0/2);
            x0_struct.S(link).S_Dy(k) = sin(theta_0/2);

            if link == 1
                if params.eta_0 == params.link(1).phi_D
                    % inner tangent between circles
                    x0_struct.S(link).alpha_P(k) = -params.eta_0*asin(2*params.link(1).R/params.L_0);
                else
                    % outer tangent between circles
                    x0_struct.S(link).alpha_P(k) = 0;
                end
            else 
                theta_prev = params.link(link-1).eta*(-omega/2 + s_k*2*omega/2);
                x0_struct.S(link).alpha_P(k) = -theta_prev/2;
            end

            x0_struct.S(link).alpha_D(k) = theta_0/2;
            x0_struct.S(link).theta(k) = theta_0;
        end
    end
end
