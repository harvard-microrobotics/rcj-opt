function [opt_vars, sol, solver] = optimize_RCJ(d_end_goal, x0_struct, params, opts)
    % build and run optimizer
    
    import casadi.*
    N_pts = params.N_pts;
    N_links = params.N_links;
    eta_0 = params.eta_0;
    L_0 = params.L_0;
    omega = params.omega;

    t = linspace(0,1,N_pts);
    theta_0 = eta_0*t*omega; % rotate_2Dation of motor
    
    % initialize
    theta_prev = theta_0; 
    eta_prev = eta_0;

    % initialize variables
    x   = {}; % decision variables
    x0 = []; % initial guess
    lbx = []; % lower decision variable bounds
    ubx = []; % upper decision variable bounds
    f   = 0; % objective
    g   = {}; % constraints
    lbg = []; % lower constraint bounds
    ubg = []; % upper constraint bounds

    for link = 1:N_links
        
        % unwrap params structure
        R = params.link(link).R; % bearing radius for this link
        L = params.link(link).L; % distance between idler bearings, or from distal frame to end-effector
        eta = params.link(link).eta; % desired rotate_2Dational direction of joint
        phi_P = params.link(link).phi_P; % proximal bearing wrap direction
        phi_D = params.link(link).phi_D; % distal bearing wrap direction
        var_bnd = params.link(link).var_bnd; % variable bounds
        const_bnd = params.link(link).const_bnd; % constraint bounds

        % storage arrays for calculated variables
        gamma_P = MX(zeros(1,N_pts));
        delta_gamma_P = MX(zeros(1,N_pts));
        c_P = MX(zeros(1,N_pts));
        gamma_D = MX(zeros(1,N_pts));
        delta_gamma_D = MX(zeros(1,N_pts));
        c_D = MX(zeros(1,N_pts));
        L_P_to_bearing = MX(zeros(1,N_pts));
        L_D_to_bearing = MX(zeros(1,N_pts));
        L_bearing = MX(zeros(1,N_pts));
        
        if link > 1
            theta_prev = vertcat(S(link-1).theta{:});
            eta_prev = params.link(link-1).eta;
        end

        for k = 1:N_pts
            % build decision variables
            S(link).S_P.x{k} = MX.sym(['S_Px_' num2str(k) '_link' num2str(link)]);
            S(link).S_P.y{k} = MX.sym(['S_Py_' num2str(k) '_link' num2str(link)]);
            S(link).S_D.x{k} = MX.sym(['S_Dx_' num2str(k) '_link' num2str(link)]);
            S(link).S_D.y{k} = MX.sym(['S_Dy_' num2str(k) '_link' num2str(link)]);
            S(link).alpha_P{k} = MX.sym(['alpha_1_' num2str(k) '_link' num2str(link)]);
            S(link).alpha_D{k} = MX.sym(['alpha_2_' num2str(k) '_link' num2str(link)]);
            S(link).theta{k} = MX.sym(['theta_'   num2str(k) '_link' num2str(link)]);
            
            % shorthand
            S_Px = S(link).S_P.x{k};
            S_Py = S(link).S_P.y{k};
            S_Dx = S(link).S_D.x{k};
            S_Dy = S(link).S_D.y{k};
            alpha_P = S(link).alpha_P{k};
            alpha_D = S(link).alpha_D{k};
            theta = S(link).theta{k};
            S_P = [S_Px; S_Py];
            S_D = [S_Dx; S_Dy];
            
            % assemble decision variables
            x = [x{:}, S_Px, S_Py, S_Dx, S_Dy, alpha_P, alpha_D, theta];

            % assemble initial guess
            x0 = [x0;
                  x0_struct.S(link).S_Px(k);
                  x0_struct.S(link).S_Py(k);
                  x0_struct.S(link).S_Dx(k);
                  x0_struct.S(link).S_Dy(k);
                  x0_struct.S(link).alpha_P(k);
                  x0_struct.S(link).alpha_D(k);
                  x0_struct.S(link).theta(k)];

            % assemble variable bounds
            lbx = [lbx; var_bnd.S_Px(1); var_bnd.S_Py(1); var_bnd.S_Dx(1); var_bnd.S_Dy(1); var_bnd.alpha_P(1); var_bnd.alpha_D(1);  var_bnd.theta(1)];
            ubx = [ubx; var_bnd.S_Px(2); var_bnd.S_Py(2); var_bnd.S_Dx(2); var_bnd.S_Dy(2); var_bnd.alpha_P(2); var_bnd.alpha_D(2);  var_bnd.theta(2)];

            % distal frame written in proximal frame
            frame_D = S_P - rotate_2D(theta)*S_D;
            S(link).frame_D{:,k} = frame_D;
            
            % transformation matrix used to calculate d_end
            T_P = [[rotate_2D(theta), frame_D+rotate_2D(theta)*[L;0]]; [0 0 1]];
            S(link).T_P{k} = T_P;

            % surface size constraints
            if params.use_constraints.S_P_size
                S_P_size = norm(S_P,2)^2;
                g   = [g(:)', {S_P_size}];
                lbg = [lbg; const_bnd.S_P_size(1)];
                ubg = [ubg; const_bnd.S_P_size(2)];
            end

            if params.use_constraints.S_D_size
                S_D_size = norm(S_D,2)^2;
                g   = [g(:)', {S_D_size}];
                lbg = [lbg; const_bnd.S_D_size(1)];
                ubg = [ubg; const_bnd.S_D_size(2)];
            end

            if k > 1
                S_P_prev = [S(link).S_P.x{k-1}; S(link).S_P.y{k-1}];
                S_D_prev = [S(link).S_D.x{k-1}; S(link).S_D.y{k-1}];
                
                % tangency constraint
                if params.use_constraints.S_tangency
                    S_P_direction = S_P - S_P_prev;
                    S_D_direction = rotate_2D(theta)*(S_D - S_D_prev);
                    S_tangency = cross2D(S_P_direction, S_D_direction);
                    g   = [g(:)', {S_tangency}];
                    lbg = [lbg; const_bnd.S_tangency(1)];
                    ubg = [ubg; const_bnd.S_tangency(2)];
                end

                % arc length constraint
                if params.use_constraints.S_arc_length
                    S_arc_length = norm(S_P - S_P_prev)^2 - norm(S_D - S_D_prev)^2;
                    g   = [g(:)', {S_arc_length}];
                    lbg = [lbg; const_bnd.S_arc_length(1)];
                    ubg = [ubg; const_bnd.S_arc_length(2)];
                end

                % surface step size constraint
                if params.use_constraints.S_step_size
                    S_step_size = norm(S_P-S_P_prev)^2;
                    g   = [g(:)', {S_step_size}];
                    lbg = [lbg; const_bnd.S_step_size(1)];
                    ubg = [ubg; const_bnd.S_step_size(2)];
                end
            end

            if k > 2
                S_P_prev_2 = [S(link).S_P.x{k-2}; S(link).S_P.y{k-2}];
                S_D_prev_2 = [S(link).S_D.x{k-2}; S(link).S_D.y{k-2}];

                % surface convexity constraints
                if params.use_constraints.S_P_convexity
                    S_P_convexity = eta * cross2D(S_P_prev - S_P_prev_2, S_P - S_P_prev_2);
                    g = [g(:)', {S_P_convexity}];
                    lbg = [lbg; const_bnd.S_P_convexity(1)];
                    ubg = [ubg; const_bnd.S_P_convexity(2)];
                end

                if params.use_constraints.S_D_convexity
                    S_D_convexity = eta * cross2D(S_D_prev - S_D_prev_2, S_D - S_D_prev_2);
                    g = [g(:)', {S_D_convexity}];
                    lbg = [lbg; const_bnd.S_D_convexity(1)];
                    ubg = [ubg; const_bnd.S_D_convexity(2)];
                end
            end

            % calculate moment arms
            if link == 1
                r_P = eta_0*(sin(alpha_P)*(L_0) + phi_P*R);
            else
                S_D_prev_link = [S(link-1).S_D.x{k}; S(link-1).S_D.y{k}];
                r_P = -eta_prev * ([-sin(alpha_P) cos(alpha_P)]*S_D_prev_link + phi_P*R);
            end
            r_D = -eta*([-sin(alpha_D) cos(alpha_D)]*S_P + phi_D*R);
            
            % moment arm constraints
            if params.use_constraints.r_P_length
                g = [g(:)', {r_P}];
                lbg = [lbg; const_bnd.r_P_length(1)];
                ubg = [ubg; const_bnd.r_P_length(2)];
            end

            if params.use_constraints.r_D_length
                g = [g(:)', {r_D}];
                lbg = [lbg; const_bnd.r_D_length(1)];
                ubg = [ubg; const_bnd.r_D_length(2)];
            end
            
            % calculate change in local angle
            if link == 1
                gamma_P(k) = alpha_P - theta_prev(k);
                c_P(k) = L_0*sin(alpha_P) + phi_P*R;
            else
                gamma_P(k) = alpha_P + theta_prev(k);
                c_P(k) = [sin(gamma_P(k)), -cos(gamma_P(k))] * S(link-1).frame_D{:,k} + phi_P*R;
            end

            gamma_D(k) = alpha_D - theta;
            c_D(k) = [-sin(alpha_D), cos(alpha_D)]*frame_D + phi_D*R;

            if k > 1
                if params.use_constraints.P_moment_arm

                    omega_1 = eta_prev*(theta_prev(k)-theta_prev(k-1));
                    omega_2 = eta*(theta-S(link).theta{k-1});

                    P_moment_arm = r_P*omega_1 + r_D*omega_2;
                    g   = [g(:)', {P_moment_arm}];
                    lbg = [lbg; const_bnd.P_moment_arm(1)];
                    ubg = [ubg; const_bnd.P_moment_arm(2)];
                end

                delta_gamma_P(k) = gamma_P(k-1) - gamma_P(k);
                delta_gamma_D(k) = gamma_D(k-1) - gamma_D(k);

                P_P_rotate_2Dation = eta_prev*delta_gamma_P(k);
                P_D_rotate_2Dation = eta*delta_gamma_D(k);
                    
                if params.use_constraints.P_P_rotate_2Dation
                    g = [g(:)', {P_P_rotate_2Dation}];
                    lbg = [lbg; const_bnd.P_P_rotate_2Dation(1)];
                    ubg = [ubg; const_bnd.P_P_rotate_2Dation(2)];
                end

                
                if params.use_constraints.P_D_rotate_2Dation
                    g = [g(:)', {P_D_rotate_2Dation}];
                    lbg = [lbg; const_bnd.P_D_rotate_2Dation(1)];
                    ubg = [ubg; const_bnd.P_D_rotate_2Dation(2)];
                end

                % Size constraints (P_P_size, P_D_size)
                if params.use_constraints.P_P_size
                    P_P_size = -(2*c_P(k-1)^2 - 4*cos(delta_gamma_P(k))*c_P(k-1)*c_P(k) + 2*c_P(k)^2)/(cos(2*delta_gamma_P(k)) - 1);
                    g   = [g(:)', {P_P_size}];
                    lbg = [lbg; const_bnd.P_P_size(1)];
                    ubg = [ubg; const_bnd.P_P_size(2)];
                end
                if params.use_constraints.P_D_size
                    P_D_size = -(2*c_D(k-1)^2 - 4*cos(delta_gamma_D(k))*c_D(k-1)*c_D(k) + 2*c_D(k)^2)/(cos(2*delta_gamma_D(k)) - 1);
                    g   = [g(:)', {P_D_size}];
                    lbg = [lbg; const_bnd.P_D_size(1)];
                    ubg = [ubg; const_bnd.P_D_size(2)];
                end
            end

            if k > 1
            % used for wire length constraint
                d_D = [-cos(alpha_D), -sin(alpha_D)] * S(link).frame_D{:,k};
                L_D_to_bearing(k) = (c_D(k-1) - c_D(k)*cos(delta_gamma_D(k)))/sin(delta_gamma_D(k)) - d_D;

                if link == 1
                    d_P = -L_0 * cos(alpha_P);
                else
                    d_P = -[cos(gamma_P(k)), sin(gamma_P(k))] * S(link-1).frame_D{:, k};
                end
                L_P_to_bearing(k) = -(c_P(k-1) - c_P(k)*cos(delta_gamma_P(k)))/sin(delta_gamma_P(k)) - d_P;

                % wire wrapped around bearing
                L_bearing(k) = R * (-phi_P*alpha_P + phi_D*alpha_D);
            end


            if k > 2
                % line intercept constraints
                P_P_intercept = (c_P(k-2) - c_P(k-1)*cos(delta_gamma_P(k-1)))/sin(delta_gamma_P(k-1)) + (c_P(k) - c_P(k-1)*cos(delta_gamma_P(k)))/sin(delta_gamma_P(k));
                if params.use_constraints.P_P_intercept
                    g = [g(:)', {P_P_intercept}];
                    lbg = [lbg; const_bnd.P_P_intercept(1)];
                    ubg = [ubg; const_bnd.P_P_intercept(2)];
                end

                P_D_intercept = ((c_D(k-2) - c_D(k-1)*cos(delta_gamma_D(k-1)))/sin(delta_gamma_D(k-1)) + (c_D(k) - c_D(k-1)*cos(delta_gamma_D(k)))/sin(delta_gamma_D(k)));
                if params.use_constraints.P_D_intercept
                    g = [g(:)', {P_D_intercept}];
                    lbg = [lbg; const_bnd.P_D_intercept(1)];
                    ubg = [ubg; const_bnd.P_D_intercept(2)];
                end

                % wire length constraint
                if params.use_constraints.L_wire
                    L_wire = L_bearing(k-1) + L_P_to_bearing(k-1) + L_D_to_bearing(k-1) + P_P_intercept - P_D_intercept - L_bearing(k) - L_P_to_bearing(k) - L_D_to_bearing(k);
                    g = [g(:)', {L_wire}];
                    lbg = [lbg; const_bnd.L_wire(1)];
                    ubg = [ubg; const_bnd.L_wire(2)];
                end

            end

        end
    end

    syms s real
    for k = 1:N_pts
        % calculate d_end using stored transformation matrices
        d_end_opt = [0; 0; 1];

        for link = wrev(1:N_links)
            d_end_opt = S(link).T_P{k}*d_end_opt;
        end
        
        % handle numeric or symbolic inputs for d_end_goal
        if isnumeric(d_end_goal)
            d_end_goal_k = d_end_goal(:,k);

        else 
            s_k = (k-1)/(N_pts-1);
            d_end_goal_k = eval(subs(d_end_goal, s, s_k));
        end

        % assemble objective
        f = f + (1/N_pts) * ((d_end_opt(1) - d_end_goal_k(1))^2 + (d_end_opt(2) - d_end_goal_k(2))^2);
    end

    % assemble and solve the nlp
    nlp = struct('f', f, 'x', vertcat(x{:}), 'g', vertcat(g{:}));
    solver = nlpsol('solver', 'ipopt', nlp, opts);
    sol = solver('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);

    % Extract optimized variables
    x_opt = full(sol.x);

    opt_vars = extract_opt_vars(x_opt, params, theta_0, d_end_goal);
end

function c = cross2D(a,b)
    % helper function to compute 2D cross product
    c = a(1)*b(2) - a(2)*b(1);
end


