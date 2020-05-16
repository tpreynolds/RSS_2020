classdef sim_data_T < handle
    
    properties
        nx(1,1) uint32
        nu(1,1) uint32
        Nsim(1,1) uint32
        tstart(1,1) double = 0
        x_nom_0(:,1) double
        x_nom_t function_handle
        u_nom_t function_handle
        t_grid(:,1) double
        Q_grid(:,:,:) double
        Y_grid(:,:,:) double
        dynamics function_handle
    end
    
    properties (SetAccess = private)
        N_sim(:,1) cell
        t_sim(:,1) cell
        x_sim(:,1) cell
        u_sim(:,1) cell
        V_sim(:,1) cell
        x0(:,1) cell
    end
    
    properties (Hidden)
        N_t(1,1) uint32 = 150
        t(1,1) double
        tspan(:,1) double
        pars(1,1) struct
    end
    
    properties (Dependent, Hidden)
       Q(:,:) double
       Y(:,:) double
       K(:,:) double
    end
    
    methods
        % basic constructor
        function obj = sim_data_T(cfnl,Nsim,tstart)
            assert(isa(cfnl,'cfga'))
            obj.nx = cfnl.nx;
            obj.nu = cfnl.nu;
            obj.Nsim = Nsim;
            if (tstart<cfnl.linear_model.t(1))
                error('tstart must be greater than %5.2f',...
                    cfnl.linear_model.t(1))
            end
            if (tstart>cfnl.linear_model.t(end))
                error('tstart must be less than %5.2f',...
                    cfnl.linear_model.t(end))
            end
            obj.tstart = tstart;
            obj.tspan  = [ tstart, cfnl.linear_model.t(end) ];
            obj.pars   = cfnl.pars;
            obj.t_grid = cfnl.linear_model.t;
            obj.Q_grid = cfnl.fnl.Q;
            obj.Y_grid = cfnl.fnl.Y;
            obj.x_nom_0 = cfnl.nominal_trj.x_t(obj.tstart);
            obj.x_nom_t = cfnl.nominal_trj.x_t;
            obj.u_nom_t = cfnl.nominal_trj.u_t;
            obj.dynamics = cfnl.dynamics;
        end
        % function to get the initial condition for a single sim
        function x_sim_0 = get_initial_condition(obj,id)
            obj.t = obj.tstart;
            Q_t = obj.Q;
            Qh  = sqrtm( Q_t );
            x0_ = zeros(obj.nx,1);
            if (id < obj.nx + 1)
                x_sim_0 = x0_ + Qh(:,id);
            elseif (id>obj.nx && id<2*obj.nx+1)
                x_sim_0 = x0_ - Qh(:,id-obj.nx);
            else
                x_sim_0 = cfga.sample_ellip(Q_t,x0_,false);
            end
            x_sim_0     = obj.x_nom_0 + x_sim_0;
            obj.x0{id}  = x_sim_0;
        end
        % function to perform on simulation and store desired data
        function one_run(obj,x_sim_0,id)
            opts = odeset('RelTol',1e-4,'AbsTol',1e-4);
            [T,X] = ode45(@(t,dx)obj.sim_dynamics(t,dx),...
                                 obj.tspan,x_sim_0,opts);
%             T = linspace(obj.tspan(1),obj.tspan(2),10000);
%             X = rk4(@(t,dx)obj.sim_dynamics(t,dx),T,x_sim_0);
            X = X(:,1:obj.nx)';
            N_T = numel(T);
            obj.N_sim{id} = N_T;
            obj.t_sim{id} = T;
            u_sim_id = zeros(obj.nu,N_T);
            V_sim_id = zeros(1,N_T);
            for k = 1:numel(T)
                obj.t         = T(k);
                dx_k          = X(:,k) - obj.x_nom_t(obj.t);
                u_sim_id(:,k) = obj.u_nom_t(obj.t) + obj.K * dx_k;
                V_sim_id(k)   = dx_k' * (obj.Q\dx_k);
            end
            obj.x_sim{id} = X;
            obj.u_sim{id} = u_sim_id;
            obj.V_sim{id} = V_sim_id;
        end
        % function to compute the differential, control, Lyap function
        function dx = sim_dynamics(obj,t,x)
            % set the time
            obj.t = t;
            x_nom = obj.x_nom_t(t);
            u_nom = obj.u_nom_t(t);
            % compute the net control
            u = u_nom + obj.K * (x - x_nom);
            % call the dynamics function
            dx = obj.dynamics(t,x,u,[],obj.pars);
        end
        % function to get the output data
        function sim_data = get_data(obj)
           sim_data = struct;
           sim_data.Nsim  = obj.Nsim;
           sim_data.tspan = [ obj.tspan(1), obj.tspan(end) ];
           sim_data.t_sim = obj.t_sim;
           sim_data.x_sim = obj.x_sim;
           sim_data.u_sim = obj.u_sim;
           sim_data.V_sim = obj.V_sim;
           sim_data.x0    = obj.x0;
        end
    end
    
     % get/set methods
    methods
        function set.t(obj,t)
            obj.t = t;
        end
        function Q = get.Q(obj)
            Q = cfga.interp_mat(obj.Q_grid,obj.t_grid,obj.t);
        end
        function Y = get.Y(obj)
            Y = cfga.interp_mat(obj.Y_grid,obj.t_grid,obj.t);
        end
        function K = get.K(obj)
            K = obj.Y/obj.Q;
        end
    end
end

