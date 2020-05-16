function [fD,dfD_dx] = pointmass_drag(x,pars)

nx = numel(x);

m = pars.m;
v = x(4:6);

speed = norm(v);
q = - 0.5 * pars.rho * pars.Sd * pars.cd;

fD = q * speed .* v;

dfD_dx = zeros(3,nx);
if (speed>1e-12)
    dfD_dx(:,4:6) = (q/m) .* ( speed * eye(3) + (v*v')/speed );
end

end

