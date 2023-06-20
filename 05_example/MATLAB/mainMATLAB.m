 % This script finds the  pareto-optimal injector inputs for pre-defined 
 % constraints.
   clear; clc; close all;

%% Initialization
 % Build and solve unconstrained ocp
   run init_ocp.m
   
 % Store unconstrained solution
   storeUnconstrainedSolution(ocp);

%% Run algorithm
 % Minimum exhaust gas temperature (Tmin)
   loadUnconstrainedSolution(ocp);
   data(1) = pareto_Tmin(ocp,parOpt,parModel,parSim);
   
 % Maximum NOx concentration (cNOx)
   loadUnconstrainedSolution(ocp);
   data(2) = pareto_cNOx(ocp,parOpt,parModel,parSim);

%% Remove temp folder with unconstrained solution
   templateFolder = what('05_example');
   tempPath = fullfile(templateFolder.path,'MATLAB','temp');
   rmpath(tempPath);
   rmdir(tempPath,'s');
 
%% Visualization
   plotResults_MATLAB(data);

   