function [f,A] = f_constraint(x,pars)

% parameters
H    = pars.H;
diff = x(pars.id_r) - pars.xc;

% value of the constraint
f = 1 - norm(H * diff); % <= 0

% partial derivative of constraint wrt state
A = zeros(1,numel(x));
A(pars.id_r) = - (diff'*(H'*H))/(norm(H*diff));

end

