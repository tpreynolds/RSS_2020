function solve_sdp(obj)

fprintf('\t\tSolving main SDP...')

M    = double(obj.M);
nx   = double(obj.nx);
nu   = double(obj.nu);
np   = double(obj.linear_model.np);
n_p  = obj.linear_model.nlg.n_p;
n_q  = obj.linear_model.nlg.n_q;
id_p = obj.linear_model.nlg.id_p;
id_q = obj.linear_model.nlg.id_q;
sze  = nx + sum(n_p) + sum(n_q);

A = obj.linear_model.A;
B = obj.linear_model.B;
C = obj.linear_model.C;
D = obj.linear_model.D;
E = obj.linear_model.E;
dt  = diff(obj.linear_model.t);
nlg = obj.linear_model.nlg.gamma;

Q_max = obj.fnl.Q_max;
R_max = obj.fnl.R_max;

Sx = obj.scale.Sx;
Su = obj.scale.Su;

lmi_tol = obj.opts.lmi_tol;
decay = obj.opts.decay_rate;
small = obj.opts.small;

% variable declarations
Q = cell(M+1,1);
Y = cell(M+1,1);
G = sdpvar(np,M,'full');    % S-procedure multipliers
rLMI = sdpvar(3,M,'full');  % relaxation parameters

% basic variable constraints
constraints = [];
for k = 1:M+1
    Q{k} = sdpvar(nx,nx);
    Y{k} = sdpvar(nu,nx);
    constraints = [ constraints, Sx*Q{k}*Sx>=0.0, ...
                                 Sx*Q{k}*Sx<=Q_max(:,:,k) ]; 
    if (k<M+1)
        for l = 1:np
           constraints = [ constraints, G(l,k)>=0.0 ];
        end
        constraints = [ constraints, rLMI(1,k)>=-1.0, rLMI(1,k)<=lmi_tol, ...
                                     rLMI(2,k)>=-1.0, rLMI(2,k)<=lmi_tol, ...
                                     rLMI(3,k)>=-1.0, rLMI(3,k)<=lmi_tol, ...
                                     rLMI(1,k)+rLMI(2,k)+rLMI(3,k)<=0.0 ];
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
    Gk  = G(:,k);
    nlg_k = nlg(:,k);
    
    R_maxk  = R_max(:,:,k);
    R_maxkp = R_max(:,:,k+1);
    
    Fkk   = Qk  * Ak'  + Ak  * Qk  + Yk'  * Bk'  + Bk  * Yk  + decay * Qk;
    Fkkp  = Qk  * Akp' + Akp * Qk  + Yk'  * Bkp' + Bkp * Yk  + decay * Qk;
    Fkpk  = Qkp * Ak'  + Ak  * Qkp + Ykp' * Bk'  + Bk  * Ykp + decay * Qkp;
    Fkpkp = Qkp * Akp' + Akp * Qkp + Ykp' * Bkp' + Bkp * Ykp + decay * Qkp;
      
    id = 1:nx;
    
    LMI1 = Fkk - dQk;
    LMI2 = Fkkp + Fkpk - 2.0 * dQk;
    LMI3 = Fkpkp - dQk;
    
    % add E terms
    for l = 1:np
        id_new = nx + id_p{l};
        Mkl    = Gk(l)*nlg_k(l)*E{l};
        % LMI1
        LMI1(id_new,id)     = Mkl';
        LMI1(id,id_new)     = Mkl;
        LMI1(id_new,id_new) = -Gk(l) * eye(n_p(l));
        % LMI2
        LMI2(id_new,id)     = 2.0*Mkl';
        LMI2(id,id_new)     = 2.0*Mkl;
        LMI2(id_new,id_new) = -2.0*Gk(l) * eye(n_p(l));
        % LMI3
        LMI3(id_new,id)     = Mkl';
        LMI3(id,id_new)     = Mkl;
        LMI3(id_new,id_new) = -Gk(l) * eye(n_p(l));
    end
    
    % add CQ+DY terms
    for l = 1:np
        id_new = nx + sum(n_p) + id_q{l};
        Nk     = C{l} * Qk  + D{l} * Yk;
        Nkp    = C{l} * Qkp + D{l} * Ykp;
        % LMI1
        LMI1(id_new,id)     = Nk;
        LMI1(id,id_new)     = Nk';
        LMI1(id_new,id_new) = -Gk(l) * eye(n_q(l));
        % LMI2
        LMI2(id_new,id)     = Nk+Nkp;
        LMI2(id,id_new)     = (Nk+Nkp)';
        LMI2(id_new,id_new) = -2.0*Gk(l) * eye(n_q(l));
        % LMI3
        LMI3(id_new,id)     = Nkp;
        LMI3(id,id_new)     = Nkp';
        LMI3(id_new,id_new) = -Gk(l) * eye(n_q(l));
    end
    
    % control constraints
    LMI4 = [ R_maxk, Yk; Yk', Qk ];
    LMI5 = [ R_maxkp, Ykp; Ykp', Qkp ];
    
    % add all to constraints
    constraints = [ constraints, LMI1 <= rLMI(1,k) * eye(sze), ...
                                 LMI2 <= rLMI(2,k) * eye(sze), ...
                                 LMI3 <= rLMI(3,k) * eye(sze), ...
                                 LMI4 >= small * eye(nx+nu), ...
                                 LMI5 >= small * eye(nx+nu) ];
    
end

% cost function
cost = -geomean(Sx*Q{1}*Sx); 

% solver options
opts = sdpsettings('solver',obj.opts.solver,'verbose',0);
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP   = 1e-5;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS     = 1e-5;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS     = 1e-5;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_MU_RED    = 1e-5;
opts.mosek.MSK_DPAR_INTPNT_CO_TOL_NEAR_REL  = 1e5;

% solve the problem
diagnose = optimize(constraints,cost,opts);

% post processing
if (diagnose.problem==0)
    obj.fnl.fnl_exist = true;
    obj.solved = true;
    for k = 1:M+1
        Qk = Sx * value(Q{k}) * Sx;
        Yk = Su * value(Y{k}) * Sx;
        % update funnel
        obj.fnl.Q(:,:,k) = Qk;
        obj.fnl.Y(:,:,k) = Yk;
    end
else
    obj.solved = false;
end

fprintf('done.\n')

end

