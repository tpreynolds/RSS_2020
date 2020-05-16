function plot_fill_ratio(obj,cfnl)

set(0,'defaultTextInterpreter',obj.interpreter,...
      'defaultLineLinewidth',obj.line_width,...
      'defaultAxesBox','on',...
      'defaultAxesFontSize',obj.font_size,...
      'defaultAxesTickLabelInterpreter',obj.interpreter,...
      'defaultLegendInterpreter',obj.interpreter,...
      'defaultLegendFontSize',obj.font_size)

pars  = cfnl.pars;
% get and scale funnel entries
Sx = diag(obj.scale_state(ones(cfnl.nx,1)));
Q     = Sx * cfnl.fnl.Q(:,:,1) * Sx;
Q_max = Sx * cfnl.fnl.Q_max(:,:,1) * Sx;

% compute fill ratio over time
Nt = 40;
t  = linspace(cfnl.linear_model.t(1),cfnl.linear_model.t(end),Nt);
fr = zeros(Nt,1);
for k = 1:Nt
    cfnl.t = t(k);
    Qt     = cfnl.Q;
    Qt_max = cfnl.Q_max;
    fr(k) = cfga.fill_ratio(Qt,Qt_max);
end

figure(3), clf
set(gcf,'Units','normalized','Position',[0.013,0.142,0.23348,0.47675])
gp = panel();
gp.pack('v',{25 []});
gp.margin = [17,15,3,3];
gp(2).pack(2,2);

% plot fill ratio
gp(1).select();
hold on, grid on, box on
plot(obj.scale_time(t),fr,'Color',obj.col.b)
set(gca,'YLim',[0,1],'XLim',obj.t_lim)
% plot position
gp(2,1,1).select()
gp(2,1,1).margin = [17,15,5,5];
hold on, grid on, box on
id     = pars.id_r;
Qr     = cfga.project_ellip(Q,id);
Qr_max = cfga.project_ellip(Q_max,id);
hr  = plot_3D_ellipsoid(Qr,zeros(3,1),obj.col.gr);
hrm = plot_3D_ellipsoid(Qr_max,zeros(3,1),obj.col.r);
title('Position [m]')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
view(3)

% plot velocity
gp(2,2,1).select()
gp(2,2,1).margin = [17,15,5,5];
hold on, grid on, box on
id     = pars.id_v;
Qv     = cfga.project_ellip(Q,id);
Qv_max = cfga.project_ellip(Q_max,id);
hv  = plot_3D_ellipsoid(Qv,zeros(3,1),obj.col.gr);
hvm = plot_3D_ellipsoid(Qv_max,zeros(3,1),obj.col.r);
title('Velocity [m/s]')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
view(3)

% plot attitude
gp(2,1,2).select()
gp(2,1,2).margin = [17,15,5,5];
hold on, grid on, box on
id     = pars.id_a;
Qa     = cfga.project_ellip(Q,id);
Qa_max = cfga.project_ellip(Q_max,id);
ha  = plot_3D_ellipsoid(Qa,zeros(3,1),obj.col.gr);
ham = plot_3D_ellipsoid(Qa_max,zeros(3,1),obj.col.r);
title('Attitude [deg]')
xlabel('roll')
ylabel('pitch')
zlabel('yaw')
view(3)

% plot angular velocity
gp(2,2,2).select()
gp(2,2,2).margin = [17,15,5,5];
hold on, grid on, box on
id     = pars.id_w;
Qw     = cfga.project_ellip(Q,id);
Qw_max = cfga.project_ellip(Q_max,id);
hw  = plot_3D_ellipsoid(Qw,zeros(3,1),obj.col.gr);
hwm = plot_3D_ellipsoid(Qw_max,zeros(3,1),obj.col.r);
title('Ang. Velocity [deg/s]')
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
view(3)

end

function [h,lims_out] = plot_3D_ellipsoid(Q,xc,col,lims_in,with2D)
    [dirs,sma] = eig(Q);
    sma = sqrt(diag(sma));
    if (det(dirs)<0)
        dirs = - dirs;
    end
    assert(det(dirs)>0)
    eul = rad2deg( dcm_2_eul(dirs) );
    
    [X,Y,Z] = ellipsoid(xc(1),xc(2),xc(3),sma(1),sma(2),sma(3),50);
    h = surf(X,Y,Z,'LineStyle','none','FaceColor',col,'FaceAlpha',0.2);
    rotate(h,[0,0,1],-eul(3));
    rotate(h,[0,1,0],-eul(2));
    rotate(h,[1,0,0],-eul(1));
      
    if (nargin>3 && ~isempty(lims_in))
        set(gca,'XLim',[-lims_in(1),lims_in(1)],...
                'YLim',[-lims_in(2),lims_in(2)],...
                'ZLim',[-lims_in(3),lims_in(3)])
    end
    
    if (nargout>1)
        lims_out = zeros(3,1);
        lims_out(1) = get(gca,'XLim');
        lims_out(2) = get(gca,'YLim');
        lims_out(3) = get(gca,'ZLim');
    end
        
    
end

