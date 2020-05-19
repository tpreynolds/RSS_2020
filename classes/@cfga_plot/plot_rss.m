function plot_rss(obj,cfnl)
%PLOT_RSS
%
%   This function creates the figures in the RSS 2020 paper

set(0,'defaultTextInterpreter',obj.interpreter,...
      'defaultLineLinewidth',obj.line_width,...
      'defaultAxesBox','on',...
      'defaultAxesFontSize',obj.font_size,...
      'defaultAxesTickLabelInterpreter',obj.interpreter,...
      'defaultLegendInterpreter',obj.interpreter,...
      'defaultLegendFontSize',obj.font_size,...
      'defaultLegendAutoUpdate','off')
col = obj.col;

nx = double(cfnl.nx);
Ix = eye(nx);

T_LIM = obj.scale_time(obj.t_lim);
X_LIM = obj.scale_state(obj.x_lim);
% manually overwrite the state limits
X_LIM(1,:) = [ -5, 5 ];
X_LIM(2,:) = [ 0, 20 ];
X_LIM(5,:) = [ -25, 25 ];
X_TIK{1} = X_LIM(1,1):2.5:X_LIM(1,2);
X_TIK{2} = X_LIM(2,1):5:X_LIM(2,2);
X_TIK{5} = X_LIM(5,1):5:X_LIM(5,2);

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

%% State funnel figure
figure(1), clf
gp = panel();
gp.margin = [19,16,3,3];
gp.pack(1,2);
gp(1,1).pack(2,1);
pos = get(gca,'Position');
wdth = 350;
hght = 350;
set(gcf,'Position',[pos(1),pos(2),wdth*2,hght]);
axis off

% make plots
for id = [1,2,5]
    if (id<3)
        gp(1,1,id,1).select();
        hold on, grid on, box on
    else
        gp(1,2).select();
        hold on, grid on, box on
    end
    
    for k = 1:N_plot-1
        tk = t_plot(k:k+1);
        Tk = T_PLOT(k:k+1);
        Xk = X_NOM(id,k:k+1);

        cfnl.t = tk(1);
        Q = cfnl.Q;
        [~,Qh_temp] = cfga.project_ellip(Q,id);
        Qh_temp = obj.scale_state(Qh_temp*Ix(:,id));
        Qh_temp = Qh_temp(id);
        cfnl.t = tk(2);
        Qp = cfnl.Q;
        [~,Qhp_temp] = cfga.project_ellip(Qp,id);
        Qhp_temp = obj.scale_state(Qhp_temp*Ix(:,id));
        Qhp_temp = Qhp_temp(id);

                    if (id==5 && k==1)
                        fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                             [ Xk(1)+Qh_temp, Xk(2)+Qhp_temp, ...
                               Xk(2)-Qhp_temp, Xk(1)-Qh_temp ],...
                             col.gr,'FaceAlpha',0.5,'LineStyle','none',...
                             'DisplayName','Funnel');
                    else
                        fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                            [ Xk(1)+Qh_temp, Xk(2)+Qhp_temp, ...
                            Xk(2)-Qhp_temp, Xk(1)-Qh_temp ],...
                            col.gr,'FaceAlpha',0.5,'LineStyle','none',...
                            'HandleVisibility','off')
                    end
    end
    % add sim data if it exists
    if (obj.sim_exist)
        for sim = 1:sim_data.Nsim
            T_SIM = sim_data.t_sim{sim};
            X_SIM = sim_data.x_sim{sim}(id,:);
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
    plot(T_PLOT,X_NOM(id,:),'Color',col.b,'DisplayName','Nominal')
    % create legend
    if (id>3)
        [~,h_leg] = legend('show','location','southwest');
        PiL = findobj(h_leg,'type','patch');
        set(PiL(1),'FaceAlpha',0.5);
    end

    set(gca,'XLim',T_LIM,'YLim',X_LIM(id,:),...
        'XTick',0:1:floor(obj.t_lim(2)),...
        'YTick',X_TIK{id})
    xlabel('Time [s]')
    ylabel(strcat('$',obj.pars.state_labels{id},'$'))
end

%% Lyapunov value figure

figure(2), clf
gp = panel();
gp.select();
hold on, grid on, box on
gp.margin = [19,16,3,3];
pos = get(gca,'Position');
set(gcf,'Position',[pos(1),pos(2),wdth,hght]);
alpha = cfnl.opts.decay_rate;
plot(T_LIM,[1,1],'--','Color',col.db,'DisplayName','Invariance')
if (obj.sim_exist)
    for sim = 1:sim_data.Nsim
        T_SIM = sim_data.t_sim{sim};
        V_SIM = sim_data.V_sim{sim};
        if (sim == 1)
            plot(T_SIM,V_SIM,'Color',col.r,'DisplayName','Test case')
            plot(T_SIM,exp(-alpha.*T_SIM),'Color',col.b,...
                    'DisplayName','Decay rate')
            legend('show','location','east')
        else
            plot(T_SIM,V_SIM,'Color',col.r,'HandleVisibility','off')
        end
    end
    set(gca,'XLim',T_LIM,'YLim',[0,1.1],...
            'XTick',0:1:floor(obj.t_lim(2)),...
            'YTick',0:0.1:1.1)
    xlabel('Time [s]')
    ylabel('$V(t)$')
end


end

