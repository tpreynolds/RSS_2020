function plot_states(obj,cfnl)

set(0,'defaultTextInterpreter',obj.interpreter,...
      'defaultLineLinewidth',obj.line_width,...
      'defaultAxesBox','on',...
      'defaultAxesFontSize',obj.font_size,...
      'defaultAxesTickLabelInterpreter',obj.interpreter,...
      'defaultLegendInterpreter',obj.interpreter,...
      'defaultLegendFontSize',obj.font_size)
col = obj.col;

nx = double(cfnl.nx);
Ix = eye(nx);
fnl_exist     = cfnl.fnl.fnl_exist;
fnl_max_exist = cfnl.fnl.fnl_max_exist;

T_LIM = obj.scale_time(obj.t_lim);
X_LIM = obj.scale_state(obj.x_lim);

N_plot = 50;
t_plot = linspace(obj.t_lim(1),obj.t_lim(2),N_plot);
T_PLOT = obj.scale_time( t_plot );

X_NOM = zeros(nx,N_plot);
for k = 1:N_plot
    X_NOM(:,k) = obj.scale_state( cfnl.nominal_trj.x_t(t_plot(k)) );
end

% if it exists, pre-scale the sim'd time and state
if (obj.sim_exist)
    sim_data = obj.sim_data;
    for sim = 1:sim_data.Nsim
        T_SIM = obj.scale_time( sim_data.t_sim{sim} );
        sim_data.t_sim{sim} = T_SIM;
        X_SIM = obj.scale_state( sim_data.x_sim{sim} );
        sim_data.x_sim{sim} = X_SIM;
    end
end

% get number of subplots
[p,~] = numSubplots(nx);

figure(1), clf
gp = panel();
gp.margin = [20,15,3,3];
gp.pack(p(1),p(2));
pos = get(gca,'Position');
wdth = 250;
hght = 250;
set(gcf,'Position',[pos(1),pos(2),wdth*p(2),hght*p(1)]);
axis off

for ii = 1:nx
    [idx,idy] = ind2sub(p,ii);
    gp(idx,idy).select();
    hold on, grid on

    for k = 1:N_plot-1
        tk = t_plot(k:k+1);
        Tk = T_PLOT(k:k+1);
        Xk = X_NOM(ii,k:k+1);
        
        if (fnl_exist)
            cfnl.t = tk(1);
            Q = cfnl.Q;
            [~,Qh_temp] = cfga.project_ellip(Q,ii);
            Qh_temp = obj.scale_state(Qh_temp*Ix(:,ii));
            Qh_temp = Qh_temp(ii);
            cfnl.t = tk(2);
            Qp = cfnl.Q;
            [~,Qhp_temp] = cfga.project_ellip(Qp,ii);
            Qhp_temp = obj.scale_state(Qhp_temp*Ix(:,ii));
            Qhp_temp = Qhp_temp(ii);
            
            if (ii==1 && k==1)
                fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                     [ Xk(1)+Qh_temp, Xk(2)+Qhp_temp, ...
                       Xk(2)-Qhp_temp, Xk(1)-Qh_temp ],...
                     col.gr,'FaceAlpha',0.5,'LineStyle','none',...
                     'DisplayName','Funnel')
            else
                fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                     [ Xk(1)+Qh_temp, Xk(2)+Qhp_temp, ...
                       Xk(2)-Qhp_temp, Xk(1)-Qh_temp ],...
                     col.gr,'FaceAlpha',0.5,'LineStyle','none',...
                     'HandleVisibility','off')
            end
             
        elseif (fnl_max_exist)
            cfnl.t = tk(1);
            Q = cfnl.Q_max;
            [~,Qh_temp] = cfga.project_ellip(Q,ii);
            Qh_temp = obj.scale_state(Qh_temp*Ix(:,ii));
            cfnl.t = tk(2);
            Qp = cfnl.Q_max;
            [~,Qhp_temp] = cfga.project_ellip(Qp,ii);
            Qhp_temp = obj.scale_state(Qhp_temp*Ix(:,ii));
            
            fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                 [ Xk(1)+Qh_temp, Xk(2)+Qhp_temp, ...
                   Xk(2)-Qhp_temp, Xk(1)-Qh_temp ],...
                 col.gr,'FaceAlpha',0.5,'LineStyle','none')
        end      

    end
    % add sim data if it exists
    if (obj.sim_exist)
        for sim = 1:sim_data.Nsim
            T_SIM = sim_data.t_sim{sim};
            X_SIM = sim_data.x_sim{sim}(ii,:);
            if (sim == 1)
                plot(T_SIM,X_SIM,'Color',col.r,...
                    'DisplayName','Test case')
            else
                plot(T_SIM,X_SIM,'Color',col.r,...
                    'HandleVisibility','off')
            end
        end
    end
    % plot nominal
    plot(T_PLOT,X_NOM(ii,:),'Color',col.b,'DisplayName','Nominal')
    if (ii==1)
        legend('show','location','southwest');
    end
    
    set(gca,'XLim',T_LIM,'YLim',X_LIM(ii,:))
    xlabel('Time [s]')
    ylabel(strcat('$',obj.pars.state_labels{ii},'$'))
    drawnow;
end
set(findall(gcf,'-property','FontSize'),'FontSize',obj.font_size)

