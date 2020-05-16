function plot_eig(obj,cfnl)

set(0,'defaultTextInterpreter',obj.interpreter,...
      'defaultLineLinewidth',obj.line_width,...
      'defaultAxesBox','on',...
      'defaultAxesFontSize',obj.font_size,...
      'defaultAxesTickLabelInterpreter',obj.interpreter,...
      'defaultLegendInterpreter',obj.interpreter,...
      'defaultLegendFontSize',obj.font_size)
col = obj.col;

N_plot = 100;
t_plot = linspace(obj.t_lim(1),obj.t_lim(2),N_plot);
T_PLOT = obj.scale_time( t_plot );
T_LIM  = obj.scale_time(obj.t_lim);

decay  = cfnl.opts.decay_rate;

Q_cond  = zeros(N_plot,1);
Q_eig   = zeros(N_plot,1);
LMI_eig = zeros(N_plot,1);

for k = 1:N_plot
   cfnl.t = t_plot(k);
   Q = cfnl.Q;
   Q_cond(k)    = cond(Q);
   Q_eig(k)     = min(eig(Q));
end

figure(4), clf
gp = panel();
gp.pack(3,1);
gp.margin = [20,20,5,7];
% plot condition number of Q(t)
gp(1,1).select();
hold on, grid on, box on
plot(T_PLOT,Q_cond,'o','Color',col.g,'MarkerFaceColor',col.g)
y_lim = get(gca,'YLim');
set(gca,'XLim',T_LIM,'YLim',[0,y_lim(2)])
xlabel('Time [s]')
ylabel('$\kappa(Q)$')

% plot min eval of Q(t)
gp(2,1).select();
hold on, grid on, box on
plot(T_PLOT,Q_eig,'Color',col.g)
plot(T_LIM,[0,0],'--','Color',[col.r,0.25])
y_lim = get(gca,'YLim');
set(gca,'XLim',T_LIM,'YLim',[0,y_lim(2)])
xlabel('Time [s]')
ylabel('Eigenvalues of Q(t)')
% plot evals of LMI(t)
gp(3,1).select();
hold on, grid on, box on
for k = 1:N_plot-1
    tk = t_plot(k:k+1);
    cfnl.t = tk(1);
    A = cfnl.A;
    B = cfnl.B;
    Y = cfnl.Y; 
    Q = cfnl.Q;
    cfnl.t  = tk(2);
    Qp = cfnl.Q;
    dQ = (1.0/diff(tk)) .* ( Qp - Q );
    LMI = A * Q + Q * A.' + B * Y + Y.' * B.' - dQ + decay * Q;
    LMI_eig(k) = max(eig(LMI));
end
plot(T_PLOT,LMI_eig,'Color',col.g)
plot(T_LIM,[0,0],'--','Color',[col.r,0.25])
set(gca,'XLim',T_LIM)
xlabel('Time [s]')
ylabel('Eigenvalues of LMI(t)')

end

