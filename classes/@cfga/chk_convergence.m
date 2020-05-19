function chk_convergence(obj)

id = { obj.pars.id_r, ...
       obj.pars.id_v, ...
	   obj.pars.id_a, ...
       obj.pars.id_w };
tol = [ obj.opts.cvrg_tol_r;
        obj.opts.cvrg_tol_v;
        obj.opts.cvrg_tol_a;
        obj.opts.cvrg_tol_w ];
nid = numel(id);

Q_max = obj.fnl.Q_max(:,:,1);
Q     = obj.fnl.Q(:,:,1);

diff = zeros(nid,1);
txt  = cell(nid,1);
cvrg = false(nid,1);

for k = 1:nid
    diff(k) = get_diff(Q_max,Q,id{k});
    if (diff(k) < tol(k))
        cvrg(k) = true;
        txt{k}  = ' ';
    else
        txt{k}  = 'x';
    end
end

% print some stuff
fprintf('\t\t r [m]  v [m/s]  a [deg]  w [deg/s]\n')
fprintf('\t\t %4.2f%s   %4.2f%s    %4.2f%s    %4.2f%s\n',...
        diff(1),txt{1},diff(2),txt{2},...
        rad2deg(diff(3)),txt{3},rad2deg(diff(4)),txt{4});

% check for convergence
if (sum(cvrg)>=obj.opts.cvrg_min)
    % did any checks not pass (ie sum(cvrg)==cvrg_min, but this is fewer
    % than nid)
    not_cvrg = find(~cvrg);
    if (any(not_cvrg))
        % if some checks did fail, check them against a slightly relaxed
        % convergence tolerance, exit if these pass since we've already
        % satisfied the minimum number of convergence checks
        temp = true;
        for k = not_cvrg
           if (diff(k) > 1.5 * tol(k))
               temp = false;
           end
        end
        obj.converged = temp;
    else
        obj.converged = true;
    end
else
    obj.converged = false;
end

% if requested, plot the fill ratio & compare funnel entries
if (obj.opts.plot_fr)
    obj.plot.plot_fill_ratio(obj)
end

end

function diff = get_diff(Q_max,Q,dims)
    Q_max_proj = cfga.project_ellip(Q_max,dims);
    Q_proj = cfga.project_ellip(Q,dims);
    Q_diff = sqrtm(Q_max_proj) - sqrtm(Q_proj);
    diff   = max(eig(Q_diff));
end
