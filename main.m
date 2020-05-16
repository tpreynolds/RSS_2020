clear; clear cfga;

run('utils/set_path.m')

cfnl = cfga('planarquad');

cfnl.plot.make_state = true;
cfnl.plot.make_control = true;
cfnl.plot.make_eig = false;
cfnl.opts.cvrg_min = 3;
cfnl.opts.max_iter = 10;
cfnl.opts.plot_fr  = false; % plot the fill ratio
cfnl.opts.decay_rate = 0.05;
cfnl.opts.lmi_tol   = 0;

cfnl.dynamics  = @dynamics;
cfnl.linearize = @linearize;

cfnl.set_scvx_data();

cfnl.synthesize_funnel(20);

cfnl.get_sim_data(12,0);

cfnl.plot.make_plots(cfnl);