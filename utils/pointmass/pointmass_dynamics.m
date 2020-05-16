function dx = pointmass_dynamics(t,x,u,ut,pars)

m = pars.m;
v = x(4:6);

% compute drag
[fD,~] = pointmass_drag(x,pars);

% interpolate control
if (~isempty(ut))
    uI = cfga.interp_vec(u,ut,t);
else
    uI = u;
end

% comput differential
dr = v;
dv = (1./m).*(uI + fD) + [ 0.0; 0.0; pars.gI];
dx = [ dr; dv ];

end
