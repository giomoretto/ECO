%% add external software packages present in the folder 04_external
% if packages are on a different path, you change the path accordingly

% acados path
acadosPath = fullfile(projectRootDir,'04_external','acados');

% casadi path
casadiVersion = 'casadi-windows-matlabR2016a-v3.4.5';
casadiPath = fullfile(projectRootDir,'04_external',casadiVersion);