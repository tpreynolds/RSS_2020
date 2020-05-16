classdef bnds_T
    
    properties
        x_min(:,1) double
        x_max(:,1) double
        u_min(:,1) double
        u_max(:,1) double
        
        dx_min(:,1) double
        dx_max(:,1) double
        du_min(:,1) double
        du_max(:,1) double
        
        u_norm_min(1,1) double
        u_norm_max(1,1) double
        
        path(1,1) struct
        
        terminal_set(:,1) double
    end
    
    properties (Hidden)
       max_tilt(1,1) double
    end
    
end

