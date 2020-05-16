classdef cfga_opts
    
    properties
        max_iter(1,1) uint32 = 10
        plot_fr(1,1) logical = false
        
        solver(1,1) string = 'mosek'
        
        decay_rate(1,1) double = 0.01
        lmi_tol(1,1) double = 1e-3
        small(1,1) double = 1e-8
        big(1,1) double = 1e4     % big enough for Yalmip, but not too big
        contract_width(1,1) double = 15
        contract_min(1,1) double = 0.6
        
        cvrg_tol_r(1,1) double = 2 % m
        cvrg_tol_v(1,1) double = 1 % m/s
        cvrg_tol_a(1,1) double = deg2rad(3) % rad
        cvrg_tol_w(1,1) double = deg2rad(2) % rad/s
        cvrg_min(1,1) double = 4
        
        nlg_method(1,1) string = 'sample'
    end
       
end

