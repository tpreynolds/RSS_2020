function set_linear_model(obj,M)

obj.M = M;

pars = obj.pars;

% initialize the linear model
obj.linear_model.init(obj);

% set the time/state/control and linearization at each temporal point
t_grid = linspace(obj.nominal_trj.t(1),obj.nominal_trj.t(end),M+1);

for k = 1:M+1
    tk      = t_grid(k);
    xk      = obj.nominal_trj.x_t(tk);
    uk      = obj.nominal_trj.u_t(tk);
    [Ak,Bk] = obj.linearize(tk,xk,uk,pars);
    
    r = rank(ctrb(Ak,Bk));
    try
        assert(r==obj.nx)
    catch
        error('Linearized model must be controllable!')
    end
    
    obj.linear_model.t(k)     = tk;
    obj.linear_model.x(:,k)   = xk;
    obj.linear_model.u(:,k)   = uk;
    obj.linear_model.A(:,:,k) = Ak;
    obj.linear_model.B(:,:,k) = Bk;
end

end