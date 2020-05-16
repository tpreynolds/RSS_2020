function contract_max_funnels(obj)

M = obj.M;
contract_min    = obj.opts.contract_min;
contract_width  = obj.opts.contract_width;
Q     = obj.fnl.Q;
Q_max = obj.fnl.Q_max;

for k = 1:M+1
    % current funnel/max entry
    Qk = Q(:,:,k);
    Q_maxk = Q_max(:,:,k);
    % compute fill ratio
    fr = cfga.fill_ratio(Qk,Q_maxk);
    % compute contraction factor
    contract = contract_min + (1.0-contract_min) ...
                    * (1.0/(1.0+exp(contract_width*(0.5-fr))));
    % contract the max funnel
    obj.fnl.Q_max(:,:,k) = contract .* Q_maxk;
end

end

