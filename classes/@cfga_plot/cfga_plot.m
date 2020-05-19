classdef cfga_plot < handle
    
    properties
        font_size(1,1) double = 16
        interpreter(1,1) string = 'latex'
        line_width(1,1) double = 1.5
        
        col(1,1) struct
        pars(1,1) struct
        
        make_rss(1,1) logical = false
        make_state(1,1) logical = false
        make_control(1,1) logical = false
        
        sim_data(1,1) struct
        sim_exist(1,1) logical = false
        
        t_lim(1,2) double
        x_lim(:,2) double
        u_lim(:,2) double
    end
    
    properties (Hidden)
       circle(2,100) double 
    end
        
    methods
        function obj = cfga_plot()
            col = struct;
            col.b  = [0,32,91]./255;
            col.g  = [10,134,61]./255;
            col.r  = [234,61,37]./255;
            col.db = [4,28,44]./255;
            col.gr = [153,153,154]./255;
            col.dr = [174,49,29]./255;
            
            col.m = [1,1,0];
            col.c = [1,0,1];
            col.cmap = [ linspace(col.m(1),col.c(1)), ...
                linspace(col.m(2),col.c(2)), ...
                linspace(col.m(3),col.c(3)) ];
            obj.col = col;
            
            angles = linspace(0,2*pi);
            obj.circle = [ cos(angles); sin(angles) ];
        end
        function init(obj,cfnl)
            obj.t_lim = [ cfnl.linear_model.t(1), cfnl.linear_model.t(end) ];
            obj.x_lim = [ cfnl.bnds.x_min, cfnl.bnds.x_max ];
            obj.u_lim = [ cfnl.bnds.u_min, cfnl.bnds.u_max ];
            obj.pars  = cfnl.pars;
        end
        function T = scale_time(~,t)
            T = t;
        end
        function X = scale_state(obj,x)
            X = x;
            [~,m] = size(X);
            for k = 1:m
                X(obj.pars.id_a,k) = rad2deg(X(obj.pars.id_a,k));
                X(obj.pars.id_w,k) = rad2deg(X(obj.pars.id_w,k));
            end
        end
        function U = scale_ctrl(~,u)
            U = u;
        end
        
        function make_plots(obj,cfnl)
           if (obj.make_state)
               obj.plot_states(cfnl)
           end
           if (obj.make_control)
               obj.plot_controls(cfnl)
           end
           if (obj.make_rss)
               obj.plot_rss(cfnl)
           end
        end
        
        plot_states(obj,cfnl);
        plot_controls(obj,cfnl);
        plot_fill_ratio(obj,cfnl);
        plot_rss(obj,cfnl);
    end
    
end

