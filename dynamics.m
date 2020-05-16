function dx = dynamics(t,x,u,ut,pars)

% parameters
m = pars.m;
J = pars.J;
l = pars.l;

% states
vI  = x(pars.id_v);
a   = x(pars.id_a);
wB  = x(pars.id_w);

if (~isempty(ut))
    uu = cfga.interp_vec(u,ut,t);
    uB = uu(pars.id_F);
else
    uB = u(pars.id_F);
end

% compute drag
[fD,~] = f_drag(x,pars);

F = [ -sin(a), -sin(a), -sin(a);
       cos(a),  cos(a),  cos(a) ];

dx = zeros(size(x));
% dynamics
dx(pars.id_r) = vI;
dx(pars.id_v) = (1./m).*( F * uB + fD ) + [ 0.0; pars.gI ];
dx(pars.id_a) = wB;
dx(pars.id_w) = (1/J).*dot([0;-l;l],uB);

end

