function synthesize_funnel(obj,M)
fprintf('CFGA solving:\n')

% set the linear model parameters
obj.set_linear_model(M);

% initialize plotting tool
obj.plot.init(obj);

% compute maximum geometric funnels
obj.get_max_funnels();

% solve the linear problem
obj.solve_sdp_no_nlg();

% main loop
for iter = 1:obj.opts.max_iter
    fprintf('\tIteration %d\n',iter)
   obj.iter = iter; 
   
   % compute nonlinear gains
   obj.get_nlg();
   
   % solve the main SDP
   obj.solve_sdp();
   
   % check for convergence
   if (obj.solved)
       obj.chk_convergence();
   end
   
   if (obj.converged)
       fprintf('converged.\n')
       break;
   else
       obj.contract_max_funnels();
   end
    
end

end

