function [A,B] = linearize(t,x,u,pars)

% sizes and parameters
nx = numel(x);
nu = numel(u);
m   = pars.m;
J   = pars.J;
l   = pars.l;

% states
a = x(pars.id_a);

% controls
uB = u(pars.id_F);

sa = sin(a);
ca = cos(a);
F = [ -sa, -sa, -sa;
       ca,  ca,  ca ];
dFda = [ -ca, -ca, -ca;
         -sa, -sa, -sa ];
[~,dfD_dx] = f_drag(x,pars);

% construct A
A = zeros(nx,nx);
A(pars.id_r,pars.id_v) = eye(2);
A(pars.id_v,pars.id_v) = dfD_dx(:,pars.id_v);
A(pars.id_v,pars.id_a) = (1/m).*dFda*uB;
A(pars.id_a,pars.id_w) = 1.0;

% construct B
B = zeros(nx,nu);
B(pars.id_v,pars.id_F) = (1/m) * F;
B(pars.id_w,pars.id_F) = (1/J) * [0, -l, l ];

end

