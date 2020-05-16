function solve_sdp_no_nlg(obj)
fprintf('\tSolving linear problem...')
% sizes and parameters
nx      = double(obj.nx);
nu      = double(obj.nu);
M       = double(obj.M);
decay   = obj.opts.decay_rate;

A       = obj.linear_model.A;
B       = obj.linear_model.B;
Q_max   = obj.fnl.Q_max;
R_max   = obj.fnl.R_max;
Sx      = obj.scale.Sx;
Su      = obj.scale.Su;
dt      = diff(obj.linear_model.t);

% solution variables
Q = cell(M+1,1);
Y = cell(M+1,1);
rLMI = sdpvar(M,3,'full');
rQ = sdpvar(M+1,1,'full');

constraints = [];
% definiteness constraints
for k = 1:M+1
    Q{k} = sdpvar(nx,nx);
    Y{k} = sdpvar(nu,nx);
    constraints = [ constraints, rQ(k)>=0.0, Sx*Q{k}*Sx>=rQ(k)*eye(nx) ];
    if (k<M+1)
       constraints = [ constraints, rLMI(k,1)>=-1.0, rLMI(k,1)<=0.0, ...
                                    rLMI(k,1)>=-1.0, rLMI(k,1)<=0.0, ...
                                    rLMI(k,1)>=-1.0, rLMI(k,1)<=0.0 ];
    end
end
% stability constraints
for k = 1:M
    Ak  = A(:,:,k);
    Akp = A(:,:,k+1);
    Bk  = B(:,:,k);
    Bkp = B(:,:,k+1);
    Qk  = Sx*Q{k}*Sx;
    Qkp = Sx*Q{k+1}*Sx;
    Yk  = Su*Y{k}*Sx;
    Ykp = Su*Y{k+1}*Sx;
    dQk = (1/dt(k)) .* ( Qkp - Qk ); 
    
    Qmaxk  = Q_max(:,:,k);
    Qmaxkp = Q_max(:,:,k+1);
    Rmaxk  = R_max(:,:,k);
    Rmaxkp = R_max(:,:,k+1);
    
    Fkk   = Qk  * Ak'  + Ak  * Qk  + Yk'  * Bk'  + Bk  * Yk  + decay * Qk;
    Fkkp  = Qk  * Akp' + Akp * Qk  + Yk'  * Bkp' + Bkp * Yk  + decay * Qk;
    Fkpk  = Qkp * Ak'  + Ak  * Qkp + Ykp' * Bk'  + Bk  * Ykp + decay * Qkp;
    Fkpkp = Qkp * Akp' + Akp * Qkp + Ykp' * Bkp' + Bkp * Ykp + decay * Qkp;
    
    LMI1 = Fkk - dQk;
    LMI2 = Fkkp + Fkpk - 2.0 * dQk;
    LMI3 = Fkpkp - dQk;
    
    LMI4 = [ Rmaxk, Yk; Yk', Qk ];
    LMI5 = [ Rmaxkp, Ykp; Ykp', Qkp ];
    
    constraints = [ constraints, LMI1 <= rLMI(k,1) * eye(nx), ...
                                 LMI2 <= rLMI(k,2) * eye(nx), ...
                                 LMI3 <= rLMI(k,3) * eye(nx), ...
                                 LMI4 >= 0.0, ...
                                 LMI5 >= 0.0, ...
                                 Qk   <= Qmaxk, ...
                                 Qkp  <= Qmaxkp ];
end

% cost function
cost = -geomean(Sx*Q{1}*Sx) - sum(rQ);
% solver options
opts = sdpsettings('solver',obj.opts.solver,'verbose',0);
% solve the problem
diagnose = optimize(constraints,cost,opts);

if (diagnose.problem==0)
    for k = 1:M+1
        Qk = cfga.cleanup_zeros( Sx * value(Q{k}) * Sx );
        Yk = cfga.cleanup_zeros( Su * value(Y{k}) * Sx );
        % compute volume ratio
        if (k==1)
            vol_ratio = det(sqrtm(Qk))/det(sqrtm(Q_max(:,:,k)));
            fprintf('volume ratio : %4.2f\n',vol_ratio)
        end
        % update max funnels
        obj.fnl.Q_max(:,:,k)  = Qk;
        obj.fnl.Qh_max(:,:,k) = sqrtm(Qk);
        obj.fnl.Y_max(:,:,k)  = Yk;
        % update linear funnels
        obj.fnl.Q_lin(:,:,k)  = Qk;
    end
end

% obj.plot.plot_states(obj);

end

