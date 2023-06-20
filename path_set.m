clear all; close all; clc;
% Set paths needed by project

% Get current folder dir
projectRootDir = fileparts(mfilename('fullpath'));


% Specify necessary folders 
includePaths = { ...
    projectRootDir; ...
    genpath(fullfile(projectRootDir, '00_utilities')); ...
    genpath(fullfile(projectRootDir, '01_model')); ...
    genpath(fullfile(projectRootDir, '02_simulation')); ...
    genpath(fullfile(projectRootDir, '03_formulation')); ...
    genpath(fullfile(projectRootDir, '05_example')); ...
    fullfile(projectRootDir, '06_build'); ...
    genpath(fullfile(projectRootDir, '06_build','acadosBuild')); ...
    fullfile(projectRootDir, '06_build','s_functions'); ...
    fullfile(projectRootDir, '06_build','s_functions','Simulation'); ...
    fullfile(projectRootDir, '06_build','s_functions','ECO'); ...
    fullfile(projectRootDir, '06_build','s_functions','ECO','c_generated_code'); ...

}; 

% Add specified folders to Matlab path
for ii = 1:size(includePaths,1)
    addpath(includePaths{ii});
end

% get path of external software packages
run path_usr.m

% include casadi folder
addpath(genpath(fullfile(casadiPath)))

% Include acados folders
addpath(fullfile(acadosPath,'interfaces','acados_matlab_octave'));
addpath(fullfile(acadosPath,'interfaces','acados_matlab_octave','acados_template_mex'));


% Set Environment Variables
env.acadosDir.name = 'ACADOS_INSTALL_DIR';
env.acadosDir.val  = fullfile(projectRootDir,'04_external','acados');
env.acadosRun.name = 'ENV_RUN';
env.acadosRun.val  = 'true';
envNames          = fieldnames(env);
for ii = 1:numel(envNames)
    setenv(env.(envNames{ii}).name,...
        env.(envNames{ii}).val);
end

% Set work directory for compiled and temporary data
workFolder = includePaths{end}; 
Simulink.fileGenControl('set', 'CacheFolder', workFolder, 'CodeGenFolder', workFolder);
save(fullfile(workFolder, 'includePaths.mat'), 'includePaths');
clear workFolder ii includePaths pathCell onPath usrPaths path_cantera path_casadi path_qpOASES

% Set default settings
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
set(groot,'DefaultTextarrowshapeInterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');

% set default colors
run RWTH_Colors.m
