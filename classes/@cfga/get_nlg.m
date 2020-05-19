function get_nlg(obj)

fprintf('\t\tComputing NLGs...')

M       = double(obj.M);
nx      = double(obj.nx);
np      = double(obj.linear_model.np);
pars    = obj.pars;
x_nom_t = obj.nominal_trj.x_t;
u_nom_t = obj.nominal_trj.u_t;

switch obj.opts.nlg_method
    
    case 'sample'
        Ns      = double(obj.linear_model.nlg.Ns);
        Ntot    = Ns * M;
        t_grid  = obj.linear_model.t;
        get_nlg = obj.linear_model.nlg.get_nlg;
        
        % loop over correction intervals
        for k = 1:M
            nlg_k = zeros(np,1);
            % sample times
            t_samp = t_grid(k) + (t_grid(k+1)-t_grid(k)).*rand(Ns);
            % loop over samples
            for samp = 1:Ns
                txt = fprintf('%02.0f%%',((Ns*(k-1)+samp)/Ntot)*100);
                obj.t = t_samp(k);
                A = obj.A;
                B = obj.B;
                Q = obj.Q_max;
                K = obj.K_max;
                x_nom = x_nom_t(obj.t);
                u_nom = u_nom_t(obj.t);
                % sample from ellipse
                dx = cfga.sample_ellip(Q,zeros(nx,1),true);
                du = K * dx;
                df = obj.dynamics(obj.t,x_nom+dx,u_nom+du,[],pars) ...
                        - obj.dynamics(obj.t,x_nom,u_nom,[],pars);
                % solve for nlg at this sample point
                [Delta,flag] = get_nlg(A,B,df,dx,du);
                % update current estimate
                if (flag==0)
                    for l = 1:np
                        if (np>1)
                            temp = norm(Delta{l},'fro');
                        else
                            temp = norm(Delta,'fro');
                        end
                        if ( temp > nlg_k(l) )
                            nlg_k(l) = temp;
                        end
                    end
                end
                fprintf(repmat('\b',1,txt));
            end
            obj.linear_model.nlg.gamma(:,k) = nlg_k;
        end
        fprintf('done.\n')
        
    case 'nlp'
        % not yet implemented
end

end

