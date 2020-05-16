function set_scvx_data(obj,filename)

fprintf('CFGA setup...')
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')

if (nargin<2 || isempty(filename))
    filename = strcat(obj.name,'_scvx');
end
       
% load data
load(strcat('data/',filename),'scvx','scales');

% problem sizes
obj.nx = scvx.nx;
obj.nu = scvx.nu;

% parameters
pars.id_r   = (1:2);
pars.id_v   = (3:4);
pars.id_a   = 5;
pars.id_w   = 6;
pars.id_F   = (1:3);
pars.gI     = scvx.auxdata.g;
pars.cd     = scvx.auxdata.cd;
pars.Sd     = scvx.auxdata.Sd;
pars.rho    = scvx.auxdata.rho;
pars.m      = scvx.auxdata.m;
pars.J      = scvx.auxdata.J;
pars.l      = scvx.auxdata.l;
pars.H      = diag([1,1]);
pars.xc     = [1;8];
pars.state_labels = {...
    'r_{\mathcal{I},x}\ [m]',...
    'r_{\mathcal{I},y}\ [m]',...
    'v_{\mathcal{I},x}\ [m/s]',...
    'v_{\mathcal{I},y}\ [m/s]',...
    '\theta\ [deg]',...
    '\omega\ [deg/s]'};
pars.control_labels = {...
    'u_{\mathcal{B},1}\ [N]',...
    'u_{\mathcal{B},2}\ [N]',...
    'u_{\mathcal{B},3}\ [N]'};
obj.pars = pars;

% set nominal trajectory
obj.nominal_trj.init(obj);
obj.nominal_trj.N = scvx.ctrl.N;
obj.nominal_trj.t = scvx.output.p .* scvx.auxdata.tau;
obj.nominal_trj.x = reshape( scvx.output.x, ...
    obj.nx, obj.nominal_trj.N );
obj.nominal_trj.u = reshape( scvx.output.u, ...
    obj.nu, obj.nominal_trj.N );
obj.nominal_trj.x_ic = obj.nominal_trj.x(:,1);
obj.nominal_trj.x_tc = obj.nominal_trj.x(:,end);
obj.nominal_trj.set_continuous_trj(pars);

% set bounds
x_min = zeros(obj.nx,1);
x_min(pars.id_r) = [ -6.0;  0.0 ];
x_min(pars.id_v) = [ -7.0; -7.0 ];
x_min(pars.id_a) = -pi/2;
x_min(pars.id_w) = -0.5;

x_max = zeros(obj.nx,1);
x_max(pars.id_r) = [ 6.0; 20.0 ];
x_max(pars.id_v) = [ 7.0;  7.0 ];
x_max(pars.id_a) = pi/2;
x_max(pars.id_w) = 0.5;

dx_min = nan(obj.nx,1);
dx_min(pars.id_a) = -deg2rad(20);

dx_max = nan(obj.nx,1);
dx_max(pars.id_a) = deg2rad(20);

u_buf = 0.1; % percent to buffer constraints
u_min = (1.0-u_buf) .* scvx.bnds.u_min;
u_max = (1.0+u_buf) .* scvx.bnds.u_max;

obj.bnds.x_min = x_min;
obj.bnds.x_max = x_max;
obj.bnds.dx_min = dx_min;
obj.bnds.dx_max = dx_max;
obj.bnds.path.x = {@(x)f_constraint(x,pars)};
obj.bnds.u_min = u_min;
obj.bnds.u_max = u_max;
obj.bnds.u_norm_min = 0;
obj.bnds.u_norm_max = norm(u_max);
obj.bnds.terminal_set = [ ...
    0.25; 0.25;
    0.15; 0.15; 
    deg2rad(3);
    deg2rad(3) ];

% scaling matrices
obj.scale.Sx = scales.Sx;
obj.scale.cx = scales.cx;
obj.scale.Su = scales.Su;
obj.scale.cu = scales.cu;

warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle')
fprintf('done.\n')
end

