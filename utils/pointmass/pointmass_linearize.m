function [A,B] = pointmass_linearize(t,x,u,pars)

nx = numel(x);
nu = numel(u);
m = pars.m;

[~,dfD_dx] = pointmass_drag(x,pars);

A           = zeros(nx,nx);
A(1:3,4:6)  = eye(3);
A(4:6,:)    = dfD_dx;

B        = zeros(nx,nu);
B(4:6,:) = (1/m)*eye(3);

end

