classdef nlg_T < handle
    
    properties
        gamma(:,:) double
        get_nlg
        n_p(:,1) double
        n_q(:,1) double
        id_p(:,1) cell
        id_q(:,1) cell
        Ns(1,1) uint32
    end
    
    methods       
        % function to initialize necessary stuff
        function init(obj,cfnl,np)
            nx = double(cfnl.nx);
            nu = double(cfnl.nu);
            np = double(np);
            C = cfnl.linear_model.C;
            D = cfnl.linear_model.D;
            E = cfnl.linear_model.E;
            % get indices and sizes
            obj.n_q = size(C{1},1);
            obj.n_p = size(E{1},2);
            obj.id_q{1} = (1:obj.n_q)';
            obj.id_p{1} = (1:obj.n_p)';
            for k = 2:np
                obj.n_q(k) = size(C{k},1);
                obj.n_p(k) = size(E{k},2);
                obj.id_q{k} = ( obj.id_q{k-1}(end)+(1:obj.n_q(k)) )';
                obj.id_p{k} = ( obj.id_p{k-1}(end)+(1:obj.n_p(k)) )';
            end
            % method specific parameters
            switch cfnl.opts.nlg_method
                case 'sample'
                    obj.Ns = 100;
                case 'nlp'
                    obj.Ns = 10;
                otherwise
                    error('unknown nonlinear gain method')
            end
            % precompile
            obj.get_nlg = obj.precompile_nlg(C,D,E,nx,nu,np);
        end
        % function to precompile the NLG computation
        function get_nlg_func = precompile_nlg(obj,C,D,E,nx,nu,np)
            % variables
            A = sdpvar(nx,nx,'full');
            B = sdpvar(nx,nu,'full');
            dx = sdpvar(nx,1,'full');
            du = sdpvar(nu,1,'full');
            df = sdpvar(nx,1,'full');
            Delta   = cell(np,1);
            cleanup = sdpvar(nx,1,'full');
            % form constraints and cost
            LHS = df - A * dx - B * du;
            RHS = zeros(nx,1);
            cost = 1e6 * norm(cleanup,2);
            for k = 1:np
                Ck = C{k};
                Dk = D{k};
                Ek = E{k};
                nqk = obj.n_q(k);
                npk = obj.n_p(k);
                Delta{k} = sdpvar(npk,nqk,'full');
                RHS = RHS + Ek * Delta{k} * ( Ck * dx + Dk * du );
                cost = cost + norm(Delta{k},'fro');
            end
            constraints = [ LHS == RHS + cleanup ]; %#ok
            % solver options
            opts = sdpsettings('solver','mosek','verbose',0);
            % parameters
            parameters  = {A,B,df,dx,du};
            % outputs
            outputs = Delta;
            % build optimizer object
            get_nlg_func = optimizer(constraints,cost,opts,...
                                                    parameters,outputs);
        end
    end
end

