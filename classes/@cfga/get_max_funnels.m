function get_max_funnels(obj)
fprintf('\tSolving for max funnels')

nx = double(obj.nx);
nu = double(obj.nu);
x_min = obj.bnds.x_min;
x_max = obj.bnds.x_max;
dx_max = obj.bnds.dx_max;
dx_min = obj.bnds.dx_min;
nx_constraints_ncvx = numel(obj.bnds.path.x);
n_bnds_x = 2*nx + sum(~isnan(dx_max)) ...
                + sum(~isnan(dx_min)) ...
                + nx_constraints_ncvx;
n_bnds_u = 2*nu;
Ix = eye(nx);

Sx = obj.scale.Sx;
Su = obj.scale.Su;

% define max state polytope
Ax = zeros(n_bnds_x,nx);
bx = zeros(n_bnds_x,1);
Ax(1:nx,:)        =  Ix;
Ax((nx+1):2*nx,:) = -Ix;
id = 2*nx + 1;
for k = 1:nx
    if (~isnan(dx_max(k)))
        Ax(id,:) = Ix(k,:);
        id = id + 1;
    end
    if (~isnan(dx_min(k)))
        Ax(id,:) = -Ix(k,:);
        id = id + 1;
    end
end
bx(1:nx,1)        =  x_max;
bx((nx+1):2*nx,1) = -x_min;

% define max control polytope
Au = [ eye(nu); -eye(nu) ];
bu = [ obj.bnds.u_max; -obj.bnds.u_min ];

% precompile optimizers before main loop
opts = sdpsettings('solver',obj.opts.solver{:},'verbose',0);

% state
xc = sdpvar(nx,1,'full');
Qh = sdpvar(nx,nx);
Qh_max = sdpvar(nx,nx);
Ax_ = sdpvar(n_bnds_x,nx,'full');
bx_ = sdpvar(n_bnds_x,1,'full');
constraints = [ Qh >= 0.0, Sx*Qh*Sx<=Qh_max ];
for ii = 1:n_bnds_x
  constraints = [ constraints, ...
      norm((Sx*Qh*Sx)*Ax_(ii,:)',2) + Ax_(ii,:)*xc <= bx_(ii) ]; %#ok
end
cost        = -geomean(Sx*Qh*Sx);
parameters  = {xc,Qh_max,Ax_,bx_};
outputs     = Sx*Qh*Sx;
get_Qh_max  = optimizer(constraints,cost,opts,parameters,outputs);

% control
uc = sdpvar(nu,1,'full');
Rh = sdpvar(nu,nu);
constraints = [ Rh>=0.0 ]; %#ok
for ii = 1:n_bnds_u
    constraints = [ constraints, ...
        norm((Su*Rh*Su)*Au(ii,:)',2) + Au(ii,:)*uc <= bu(ii) ]; %#ok
end
cost       = -geomean(Su*Rh*Su);
parameters = uc;
outputs    = Su*Rh*Su;
get_Rh_max = optimizer(constraints,cost,opts,parameters,outputs);
       
% loop over correction nodes
for k = 1:obj.M+1
    
    % solve for max state funnel
    xc = obj.linear_model.x(:,k);
    if (k==obj.M+1)
        Qh_max = diag(obj.bnds.terminal_set);
    else
        Qh_max = obj.opts.big * eye(nx); 
    end
    id = 2*nx+1;
    for ii = 1:nx
        if (~isnan(dx_max(ii)))
            bx(id,1) = xc(ii) + dx_max(ii);
            id = id + 1;
        end
        if (~isnan(dx_min(ii)))
            bx(id,1) = -(xc(ii) + dx_min(ii));
            id = id + 1;
        end
    end
    for cnstr = 1:nx_constraints_ncvx
        [fcnstr,Acnstr] = obj.bnds.path.x{cnstr}(xc);
        Ax(id,:) = Acnstr;
        bx(id)   = Acnstr*xc - fcnstr;
        id = id + 1;
    end
    [Qh,flag] = get_Qh_max(xc,Qh_max,Ax,bx);
    if ~(flag==0)
        error('max state funnel infeasible at iteration %d',k)
    else
        Qh = cfga.cleanup_zeros(Qh);
        obj.fnl.Qh_max(:,:,k) = Qh;
        obj.fnl.Q_max(:,:,k)  = Qh * Qh;
    end
    
    % solve for max control funnel
    uc = obj.linear_model.u(:,k);
    [Rh,flag] = get_Rh_max(uc);
    if ~(flag==0)
        error('max control funnel infeasible at iteration %d',k)
    else
        Rh = cfga.cleanup_zeros(Rh);
        obj.fnl.Rh_max(:,:,k) = Rh;
        obj.fnl.R_max(:,:,k)  = Rh * Rh;
    end
    
    fprintf('.')
end
% set existence of max funnel to true
obj.fnl.fnl_max_exist = true;

fprintf('done.\n')

end