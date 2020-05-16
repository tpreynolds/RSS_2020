function plot_controls(obj,cfnl)

set(0,'defaultTextInterpreter',obj.interpreter,...
      'defaultLineLinewidth',obj.line_width,...
      'defaultAxesBox','on',...
      'defaultAxesFontSize',obj.font_size,...
      'defaultAxesTickLabelInterpreter',obj.interpreter,...
      'defaultLegendInterpreter',obj.interpreter,...
      'defaultLegendFontSize',obj.font_size)
col = obj.col;

nu = double(cfnl.nu);
Iu = eye(nu);
fnl_max_exist = cfnl.fnl.fnl_max_exist;

T_LIM = obj.scale_time(obj.t_lim);
U_LIM = obj.scale_ctrl(obj.u_lim);

N_plot = 50;
t_plot = linspace(obj.t_lim(1),obj.t_lim(2),N_plot);
T_PLOT = obj.scale_time( t_plot );

U_NOM = zeros(nu,N_plot);
for k = 1:N_plot
    U_NOM(:,k) = obj.scale_ctrl( cfnl.nominal_trj.u_t(t_plot(k)) );
end

% if it exists, pre-scale the sim'd time and control
if (obj.sim_exist)
    sim_data = obj.sim_data;
    for sim = 1:sim_data.Nsim
        T_SIM = obj.scale_time( sim_data.t_sim{sim} );
        sim_data.t_sim{sim} = T_SIM;
        U_SIM = obj.scale_ctrl( sim_data.u_sim{sim} );
        sim_data.u_sim{sim} = U_SIM;
    end
end

[p,~] = numSubplots(nu);

figure(2), clf
gp = panel();
gp.margin = [25,17,5,7];
gp.pack(p(1),p(2));
pos = get(gca,'Position');
wdth = 300;
hght = 300;
set(gcf,'Position',[pos(1)+250,pos(2)+250,wdth*p(2),hght*p(1)]);
axis off

for ctrl = 1:nu
    [ii,jj] = ind2sub(p,ctrl);
    gp(ii,jj).select();
    hold on, grid on, box on
    
    for k = 1:N_plot-1
        tk = t_plot(k:k+1);
        Tk = T_PLOT(k:k+1);
        Uk = U_NOM(ctrl,k:k+1);
        
        if (fnl_max_exist)
            cfnl.t = tk(1);
            R  = cfnl.R_max;
            [~,Rh_temp] = cfga.project_ellip(R,ctrl);
            Rh_temp     = obj.scale_ctrl(Rh_temp*Iu(:,ctrl));
            cfnl.t = tk(2);
            Rp = cfnl.R_max;
            [~,Rhp_temp] = cfga.project_ellip(Rp,ctrl);
            Rhp_temp = obj.scale_ctrl(Rhp_temp*Iu(:,ctrl));
            
            fill([ Tk(1), Tk(2), Tk(2), Tk(1) ],...
                 [ Uk(1)+Rh_temp,  Uk(2)+Rhp_temp, ...
                   Uk(2)-Rhp_temp, Uk(1)-Rh_temp ],...
                 col.gr,'FaceAlpha',0.5,'LineStyle','none',...
                 'HandleVisibility','off')
        end
    end
    % add sim data if it exists
    if (obj.sim_exist)
        for sim = 1:sim_data.Nsim
            T_SIM = sim_data.t_sim{sim};
            U_SIM = sim_data.u_sim{sim}(ctrl,:);
            if (sim == 1)
                plot(T_SIM,U_SIM,'Color',col.r,...
                    'DisplayName','Test case')
            else
                plot(T_SIM,U_SIM,'Color',col.r,...
                    'HandleVisibility','off')
            end
        end
    end
    % plot nominal
    plot(T_PLOT,U_NOM(ctrl,:),'Color',col.b,'DisplayName','Nominal')
    
    if (ctrl==1)
        legend('show','location','best')
    end
    
    set(gca,'XLim',T_LIM,'YLim',U_LIM(ctrl,:))
    xlabel('Time [s]')
    ylabel(strcat('$',obj.pars.control_labels{ctrl},'$'))
    drawnow;
end


end

