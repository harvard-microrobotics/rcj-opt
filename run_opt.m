% add casadi to path and set enviornment before use

clear
close all

folderPath = ''; % path to save animation mp4

N_pts = 30; % number of time steps
N_links = 2; % number of links

% generate parameters and bounds
params = generate_params(N_links, N_pts);

% generate default goal path
dend_goal = generate_goal_path(params);

% generate initial guess
x0_struct = generate_x0(params);

% set optimization parameters
opts = struct;
opts.ipopt.max_iter = 500;

% build and run optimization
[opt_vars_noDend, sol_noDend, solver_noDend] = optimize_RCJ(dend_goal, x0_struct, params, opts);

% visualize results
filename = 'RCJ_animation';
animate_RCJ(opt_vars_noDend, params, [folderPath, filename]);

