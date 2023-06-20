% This script runs the solver to find optimal injector inputs for a
% pre-defined reference trajectory. You can either use the template or
% define your own trajectory.
clear; clc; close all;
load parECO.mat

%% Simulation Settings
 % Measurement noise
   parECO.enNoise	 = true;
 % Disturbances
   parECO.enDist	 = true;
 % FB Controllers
   parECO.enFB.imep  = true;
   parECO.enFB.dpMax = true;
   parECO.enFB.Tevo  = true;
   parECO.enFB.cNOx  = true;

%% Trajectory definition
 % Template
   parECO.ref.dT	 = 5; % new reference after dT seconds
   parECO.ref.imep	 = [  6   6   9   9   7   7 ] * 1e5;	% [Pa]
   parECO.ref.dpMax	 = [ 15   4   4   4   4   4 ] * 1e5;	% [Pa]
   parECO.ref.Tevo   = [  0   0   0 680   0   0 ] + 273;	% [K]
   parECO.ref.cNOx	 = [100 100 100 100 100  11 ] * 1e2;	% [ppm]
 % Custom
   % parECO.ref.dT     = ...;
   % parECO.ref.imep   = [ ... ] * 1e5;
   % parECO.ref.dpMax  = [ ... ] * 1e5;
   % parECO.ref.Tevo   = [ ... ] + 273;
   % parECO.ref.cNOx   = [ ... ] * 1e2;
   
 % The algorithm is sensitive towards initial conditions. To ensure correct
 % behaviour of the algorithm, initial conditions are fixed
   parECO.ref.time  = 0 : parECO.ref.dT : parECO.ref.dT*(numel(parECO.ref.imep)-1);
   parECO.ref.imep  = [parECO.Reference.IMEP	parECO.ref.imep(1)  parECO.ref.imep ];
   parECO.ref.dpMax = [parECO.Reference.dpMax	parECO.ref.dpMax(1) parECO.ref.dpMax];
   parECO.ref.Tevo  = [parECO.Reference.Tmin	parECO.ref.Tevo(1)  parECO.ref.Tevo];
   parECO.ref.cNOx  = [parECO.Reference.cNOx	parECO.ref.cNOx(1)  parECO.ref.cNOx];
   parECO.ref.time  = [0 1 parECO.ref.time + parECO.ref.dT];
   
%% Run Simulink simulation
   fprintf('Run Simulink Simulation.. ');
   simOut    = sim('ECO_openSource_2017','SaveOutput','on','SrcWorkspace','current');
   fprintf('Done\n\n');
   
%% Visualization
   plotResults_Simulink(parECO.ref, simOut, parECO);
   