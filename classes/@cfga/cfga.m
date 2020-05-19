classdef cfga < handle
    
    properties (Access = public)
        name(1,1) string
        pars(1,1) struct
        opts(1,1) cfga_opts
        plot(1,1) cfga_plot
        dynamics function_handle
        linearize function_handle
    end
    
    properties (SetAccess = private)
        nx(1,1) uint32
        nu(1,1) uint32
        M(1,1) uint32
        fnl(1,1) fnl_T
        linear_model(1,1) linear_model_T
        nominal_trj(1,1) nominal_trj_T
        bnds(1,1) bnds_T
    end
    
    properties (Access = public, Hidden)
        t(1,1) double
    end
    
    properties (SetAccess = private, Hidden)
        iter(1,1) uint32
        scale(1,1) struct
        converged(1,1) logical
        solved(1,1) logical
    end
    
    properties (Dependent, Hidden)
        A(:,:) double
        B(:,:) double
        Q(:,:) double
        Y(:,:) double
        Q_max(:,:) double
        Q_lin(:,:) double
        Y_max(:,:) double
        R_max(:,:) double
        K(:,:) double
        K_max(:,:) double
    end
    
    methods (Access = public)
        % basic constructor
        function obj = cfga(name)
            if (isa(name,'char'))
                obj.name = name;
            else
                error('Invalid input, must be a string')
            end
        end
        % function to set scvx data from file
        set_scvx_data(obj,filename)
        % function to set linear model data
        set_linear_model(obj,M)
        % function to synthesize funnel
        synthesize_funnel(obj,M)
        % function to get maximum geometric funnels
        get_max_funnels(obj)
        % function to solve the linear-SDP problem
        solve_sdp_no_nlg(obj)
        % function to compute nonlinear gains
        get_nlg(obj)
        % function to solve the main SDP
        solve_sdp(obj)
        % function to check convergence
        chk_convergence(obj)
        % function to contract max funnels
        contract_max_funnels(obj)
        % function to simulate various initial conditions from funnel
        get_sim_data(obj,Nsim,tstart)
        % function to save the minimum data to construct any trj
        save_data(obj,save_all)
        % function to print out the funnel size at time t
        function print_funnel_size(obj,t)
            obj.t = t;
            Qt = obj.Q;
            fprintf('Radii of funnel projections onto each state dimension')
            fprintf('at time t=%4.2f:\n',t)
            Ix = eye(obj.nx);
            for dim = 1:obj.nx
                [~,Qdim] = cfga.project_ellip(Qt,dim);
                Qdim = obj.plot.scale_state(Qdim*Ix(:,dim));
                dim_label = obj.pars.state_labels{dim}; 
                fprintf('%s : %5.2f\n',dim_label,Qdim(dim));
            end
        end
    end
       
    % static methods
    methods (Static)
        % function to map a vector to 3x3 cross product matrix
        function v_x = skew3(v)
            if (numel(v)~=3)
                error('input vector must have 3 elements')
            end
            v_x = [ 0.0, -v(3), v(2);
                    v(3), 0.0, -v(1);
                   -v(2), v(1), 0.0 ]; 
        end
        % function to get current temporal interval
        function interval = get_interval(t_grid,t)
            if (t<t_grid(1))
                error('time %f is less than the grid lower value of %f',...
                    t,t_grid(1));
            end
            if (t>t_grid(end))
                error('time %f is greater than the grid upper value of %f',...
                    t,t_grid(end));
            end
            interval = sum(t>t_grid);
            if (interval==0)
                interval = 1;
            end
        end
        % function to linearly interpolate vectors
        function v_t = interp_vec(v,t_grid,t)
            if (numel(t_grid)~=size(v,2))
                error('size of input vector must match number of temporal points')
            end
            k = cfga.get_interval(t_grid,t);
            sigma = (t_grid(k+1)-t)/(t_grid(k+1)-t_grid(k));
            v_t = sigma.*v(:,k) + (1-sigma).*v(:,k+1);
        end
        % function to linearly interpolate matrices
        function M_t = interp_mat(M,t_grid,t)
            if (numel(t_grid)~=size(M,3))
                error('size of input vector must match number of temporal points')
            end
            k = cfga.get_interval(t_grid,t);
            sigma = (t_grid(k+1)-t)/(t_grid(k+1)-t_grid(k));
            M_t = sigma.*M(:,:,k) + (1-sigma).*M(:,:,k+1);
        end
        % function to project an ellipse onto specified dimensions
        function [Q_proj,Qh_proj] = project_ellip(Q,dims)
            n = size(Q,1);
            ndims = numel(dims);
            if (ndims>n)
                error('requested dimensions not compatible with ellipse')
            end
            % build projection matrix
            In = eye(n);
            Im = eye(ndims);
            T  = zeros(n,ndims);
            for k = 1:ndims
                T(:,k) = In(:,dims(k));
            end
            P = In/Q;
            % the ldl decomposition here is more robust than the cholesky
            % since Yalmip can produce some solutions that have a very
            % small negative eigenvalues. Mathematically this procedure is
            % identical, but ldl does not throw an error if P is not pos
            % def. Moreover, the factor L is not necessarily lower
            % triangular since ldl does some pivoting that a cholesky
            % factorization does not need to do. However, the results are
            % the same.
            % L = chol(P,'lower');
            [L_,D_,P_] = ldl(P);
            L = P_ * L_ * sqrtm(D_);
            A = T' * (In/L)';
            [U,S,~] = svd(A,'econ');
            iS = Im/S;
            P_proj  = U * iS * iS * U';
            Q_proj  = Im/P_proj;
            Qh_proj = sqrtm(Q_proj);
        end
        % function to sample uniformly from an ellipse
        function xs = sample_ellip(Q,xc,bdd)
            try assert(isa(bdd,'logical')); catch
                error('third input must be a ''logical'' type')
            end
            n = size(Q,1);
            if (numel(xc)~=n)
                error('dimensions of inputs Q and xc must match')
            end
            switch bdd
                case true
                    % uniformly sample unit sphere 
                    z = randn(n,1);
                    z = z./norm(z);
                    % map to ellipse
                    xs = xc + sqrtm(Q) * z;
                case false
                    % compute bounding box
                    bbox = zeros(n,2);
                    for k = 1:n
                        [~,bbox(k,2)] = cfga.project_ellip(Q,k);
                    end
                    bbox(:,1) = - bbox(:,2);
                    outside = true;
                    while outside 
                        % sample
                        xs = bbox(:,1) + (bbox(:,2)-bbox(:,1)).*rand(n,1);
                        % test
                        val = (xs-xc)'*(Q\(xs-xc));
                        if (val<1.001)
                            outside = false;
                        end
                    end
            end
            
            
        end
        % function to compute the fill ratio between two ellipsoids
        function val = fill_ratio(Q,Q_max)
            n = size(Q,1);
            val = 1e6;
            for k = 1:n
               [~,Qh_proj_k] = cfga.project_ellip(Q,k);
               [~,Qh_max_proj_k] = cfga.project_ellip(Q_max,k);
               temp = Qh_proj_k/Qh_max_proj_k;
               if (temp<val)
                   val = temp;
               end
            end
        end
        % function to clean up annoyingly small numbers 
        function M = cleanup_zeros(M)
           [n,m] = size(M);
           for row = 1:n
               for col = 1:m
                   if (abs(M(row,col))<eps)
                       M(row,col) = 0.0;
                   end
               end
           end
        end
    end

    % get/set methods
    methods
        function set.t(obj,t)
            if (t<0)
                error('Time must be nonnegative')
            end
            if (t>obj.linear_model.t(end))
                error('Time must be < the final time')
            end
            obj.t = t;
        end
        function A = get.A(obj)
            A_grid = obj.linear_model.A;
            t_grid = obj.linear_model.t;
            A = cfga.interp_mat(A_grid,t_grid,obj.t);
        end
        function B = get.B(obj)
            B_grid = obj.linear_model.B;
            t_grid = obj.linear_model.t;
            B = cfga.interp_mat(B_grid,t_grid,obj.t);
        end
        function Q = get.Q(obj)
            Q_grid = obj.fnl.Q;
            t_grid = obj.linear_model.t;
            Q = cfga.interp_mat(Q_grid,t_grid,obj.t);
        end
        function Y = get.Y(obj)
            Y_grid = obj.fnl.Y;
            t_grid = obj.linear_model.t;
            Y = cfga.interp_mat(Y_grid,t_grid,obj.t);
        end
        function K = get.K(obj)
            K = obj.Y/obj.Q;
        end
        function Q_max = get.Q_max(obj)
            Q_max_grid = obj.fnl.Q_max;
            t_grid = obj.linear_model.t;
            Q_max  = cfga.interp_mat(Q_max_grid,t_grid,obj.t);
        end
        function Q_lin = get.Q_lin(obj)
            Q_lin_grid = obj.fnl.Q_lin;
            t_grid = obj.linear_model.t;
            Q_lin  = cfga.interp_mat(Q_lin_grid,t_grid,obj.t);
        end
        function Y_max = get.Y_max(obj)
            Y_max_grid = obj.fnl.Y_max;
            t_grid = obj.linear_model.t;
            Y_max  = cfga.interp_mat(Y_max_grid,t_grid,obj.t);
        end
        function R_max = get.R_max(obj)
           R_max_grid = obj.fnl.R_max;
           t_grid = obj.linear_model.t;
           R_max = cfga.interp_mat(R_max_grid,t_grid,obj.t);
        end
        function K_lin = get.K_max(obj)
           K_lin = obj.Y_max/obj.Q_lin; 
        end
    end
end

