clear; clear cfga;

run('utils/set_path.m')

% instantiate a CFGA object
cfnl = cfga('planarquad');

% general options
cfnl.opts.cvrg_min      = 3;
cfnl.opts.max_iter      = 10;
cfnl.opts.decay_rate    = 0.01;
cfnl.opts.lmi_tol       = 0;
cfnl.plot.make_rss      = true;

% variable contraction parameters
cfnl.opts.contract_min      = 0.7;
cfnl.opts.contract_width    = 15;

% convergence tolerances
cfnl.opts.cvrg_tol_a = deg2rad(3);
cfnl.opts.cvrg_tol_w = deg2rad(3);

% attach dynamics and linearization functions
cfnl.dynamics  = @dynamics;
cfnl.linearize = @linearize;

% set SCvx data from stored trajectory
cfnl.set_scvx_data();

% synthesize funnels
cfnl.synthesize_funnel(20);

% simulate some test cases ( # cases, t_{start} )
cfnl.get_sim_data(25,0);

% make the desired plots
cfnl.plot.make_plots(cfnl);