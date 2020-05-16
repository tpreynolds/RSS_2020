function save_data(obj,save_all)

filename_full = strcat('data/',obj.name,'_full');
filename_min = strcat('data/',obj.name,'_min');

Q_grid  = obj.fnl.Q;
Y_grid  = obj.fnl.Y;
u_nom   = obj.nominal_trj.u;
x_nom_0 = obj.nominal_trj.x(:,1);

% save all data if requested
if (save_all)
    save(filename_full{:},'obj')
end

% save minimum data always
save(filename_min{:},'Q_grid','Y_grid','u_nom','x_nom_0');

% get number of bytes required for minimum data
nB = 0;
nB = add_var_size(nB,Q_grid);
nB = add_var_size(nB,Y_grid);
nB = add_var_size(nB,u_nom);
nB = add_var_size(nB,x_nom_0);

fprintf('CFGA requires %5.2f kB of data\n',nB*1e-3)

end

function nB = add_var_size(nB,var)
    switch class(var)
        case 'double'
            B_per_el = 8;
        case 'single'
            B_per_el = 4;
        case {'int8','uint8'}
            B_per_el = 1;
        case {'int16','uint16'}
            B_per_el = 2;
        case {'int32','uint32'}
            B_per_el = 4;
        case {'int64','uint64'}
            B_per_el = 8;
        case 'logical'
            B_per_el = 1;
        case 'char'
            B_per_el = 2;
        otherwise
            error('unrecognized data type')
    end
    
    nB = nB + numel(var) * B_per_el;
end

