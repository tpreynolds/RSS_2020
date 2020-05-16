classdef nominal_trj_T < handle
    
    properties
        N(1,1) uint32
        t(:,:) double
        x(:,:) double
        u(:,:) double
        x_t function_handle
        u_t function_handle
    end
    
    properties (Hidden)
       x_ic(:,1) double
       x_tc(:,1) double
       dynamics function_handle
    end
    
    methods
        
        function init(obj,cfnl)
            obj.dynamics = cfnl.dynamics;
        end
        
        function set_continuous_trj(obj,pars)
            t_grid  = obj.t;
            tspan   = [obj.t(1), obj.t(end)];
            opts    = odeset('RelTol',1e-6,'AbsTol',1e-6);
            x0      = obj.x_ic;
            [T,X]   = ode45(@(t,x)obj.dynamics(t,x,obj.u,t_grid,pars),...
                        tspan,x0,opts);
            obj.x_t = @(t)cfga.interp_vec(X',T,t);
            obj.u_t = @(t)cfga.interp_vec(obj.u,t_grid,t);
        end
    end
    
end

