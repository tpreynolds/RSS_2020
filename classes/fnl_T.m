classdef fnl_T
    
    properties
        Q(:,:,:) double
        Qh(:,:,:) double
        Y(:,:,:) double
        Q_max(:,:,:) double
        Qh_max(:,:,:) double
        R_max(:,:,:) double
        Rh_max(:,:,:) double
        Y_max(:,:,:) double
        Q_lin(:,:,:) double
        Y_lin(:,:,:) double
        fnl_exist(1,1) logical = false
        fnl_max_exist(1,1) logical = false
    end
    
end

