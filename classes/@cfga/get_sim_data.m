function get_sim_data(obj,Nsim,tstart)

fprintf('CFGA: generating simulation data...')

% instantiate a sim object
sim_obj = sim_data_T(obj,Nsim,tstart);

for sm = 1:Nsim
    txt = fprintf('%02.0f%%',(sm/Nsim)*100);
    
    % get initial condition
    x0 = sim_obj.get_initial_condition(sm);
    
    % simulate 
    sim_obj.one_run(x0,sm);
        
    fprintf(repmat('\b',1,txt))
end

% extract necessary data and let plotting tools know that sim_data exists
obj.plot.sim_data  = sim_obj.get_data;
obj.plot.sim_exist = true;

fprintf('done\n')
end