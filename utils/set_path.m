path(pathdef)

% add the CFGA files
addpath(genpath('../utils/'))
addpath('../classes/')

% add the Yalmip files
try
    yalmiproot = 'YALMIP/';
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

% add the solver files
try
    solverroot = '../../Desktop/mosek/9.1/toolbox/r2015a/';
    addpath(genpath(solverroot))
    clear solverroot
catch
    fprintf('Did not find a SOLVER installation in the specified')
    fprintf(' directory, stopping')
end