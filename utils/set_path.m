path(pathdef)

% add the CFGA files
addpath(genpath('../utils/'))
addpath('../classes/')

% add the Yalmip files
try
    yalmiproot = '../../Desktop/YALMIP-master/';
    addpath(yalmiproot)
    addpath(strcat(yalmiproot,'extras'))
    addpath(strcat(yalmiproot,'solvers'))
    addpath(strcat(yalmiproot,'modules'))
    addpath(strcat(yalmiproot,'modules/parametric'))
    addpath(strcat(yalmiproot,'modules/moment'))
    addpath(strcat(yalmiproot,'modules/global'))
    addpath(strcat(yalmiproot,'modules/sos'))
    addpath(strcat(yalmiproot,'operators'))
    clear yalmiproot 
catch 
    fprintf('Did not find a YALMIP installation on the Desktop, stopping')
    return;
end

% add the mosek files
try
    addpath(genpath('../../Desktop/mosek/9.1/toolbox/r2015a/'))
catch
    fprintf('Did not find a MOSEK installation on the Desktop, stopping')
end
