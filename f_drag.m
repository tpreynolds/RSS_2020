function [fD,dfD_dx] = f_drag(x,pars)

% paramters
m = pars.m;
q = - 0.5 * pars.rho * pars.Sd * pars.cd;

% states
v = x(pars.id_v);
speed = norm(v);

% compute drag force
fD = q * speed .* v;

% compute partial of drag with respect to the state
nv = numel(pars.id_v);
dfD_dx = zeros(nv,numel(x));
if (speed>1e-12)
    dfD_dx(:,pars.id_v) = (q/m) .* ( speed * eye(nv) + (v*v')/speed );
end

end

