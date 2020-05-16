classdef linear_model_T < handle
    
    properties
        A(:,:,:) double
        B(:,:,:) double
        C(:,1) cell
        D(:,1) cell
        E(:,1) cell
        t(1,:) double
        x(:,:) double
        u(:,:) double
        np(1,1) uint32
        nlg(1,1) nlg_T
    end
    
    methods 
        function init(obj,cfnl)
            assert(isa(cfnl,'cfga'))
            % sizes and parameters
            nx = cfnl.nx;
            nu = cfnl.nu;
            id_v = cfnl.pars.id_v;
            id_a = cfnl.pars.id_a;
            id_w = cfnl.pars.id_w;
            
            % velocity to velocity
            obj.add_state_nonlinearity(id_v,id_v,nx,nu);
            
            % attitude to velocity
            obj.add_state_nonlinearity(id_a,id_v,nx,nu);
            
%             % control to velocity
%             obj.add_control_nonlinearity(id_v,nx,nu);
%             
%             % control to angular velocity
%             obj.add_control_nonlinearity(id_w,nx,nu);
%             
            % add the number of nonlinear channels
            obj.np = numel(obj.C);
            
            % initialize nonlinear gain
            obj.nlg.init(cfnl,obj.np);
        end
                
        function add_state_nonlinearity(obj,id_from,id_to,nx,nu)
            % set the sizes
            C_ = zeros(numel(id_from),nx);
            D_ = zeros(numel(id_from),nu);
            E_ = zeros(nx,numel(id_to));
            % set non-zero blocks
            C_(:,id_from) = eye(numel(id_from));
            E_(id_to,:)   = ones(numel(id_to));
            % append to current object
            obj.C{end+1,1} = C_;
            obj.D{end+1,1} = D_;
            obj.E{end+1,1} = E_;
        end
        
        function add_control_nonlinearity(obj,id_to,nx,nu)
            % set the sizes
            C_ = zeros(nu,nx);
            E_ = zeros(nx,numel(id_to));
            % set non-zero blocks
            D_          = eye(nu);
            E_(id_to,:) = ones(numel(id_to));
            % append to current object
            obj.C{end+1,1} = C_;
            obj.D{end+1,1} = D_;
            obj.E{end+1,1} = E_;
        end
    end
    
end

