function [] = getSFun(sqpPrepFun)
%% -----------------------------------------------------------------------%
%       Create CasADi functions for NLP with SQP & Multiple Shooting      %
%-------------------------------------------------------------------------%

%% Preliminary Code
% Get build folder path
  buildFolder = what('06_build');
  sFunFolder  = fullfile(buildFolder.path,'s_functions');

% If there exists a folder already, remove from path and delete it
  pathCell = regexp(path, pathsep, 'split');
  if ispc
      onPath = any(strcmpi(sFunFolder, pathCell));
  else
      onPath = any(strcmp(sFunFolder, pathCell));
  end

%% Run S-Function Generation
% Compile S-functions, needed for the Simulink simulation
fprintf('\nGenerating S-function for Simulation... \n\n');
% The following function creates the mex files needed by the S-functions
% in the simulink environment and saves them in the folder "sFunctions".
% Furhter, for every mex file, C code is generated. This code will be
% later on used to run the derived simulation on a real-time-processor.
%
% This function works with the MinGW C (not C++) Compiler to generate mex
% files. It can be installed under HOME-->Add-Ons. After installation the
% mex setup has to be configured accordingly. This can be done in the
% command window with the command "mex -setup".

% If SFun Folder is on path, remove it
  if onPath
      rmpath(sFunFolder);
  end
% Clear SFun Folder
  clear mex; %#ok<CLMEX>
  if isfolder(sFunFolder)
      rmdir(sFunFolder,'s')
  end
  pause(0.2);

% For all functions, run generation
  fieldNames = fieldnames(sqpPrepFun);
  for i = 1:numel(fieldNames)
      % Define codegeneration options
      OptionsCg.output_is_sparse     = 1;
      OptionsCg.build_directory      = fullfile(sFunFolder,'Simulation');
      OptionsCg.output_is_matrix     = 0;
      OptionsCg.common_sourcefile    = fieldNames{i};
      OptionsCg.generate_rtwmakecfg  = 1;
      OptionsCg.enable_tictoc_dspace = 0;
      OptionsCg.options_s_function   = {'SS_OPTION_CAN_BE_CALLED_CONDITIONALLY', ...
                                        'SS_OPTION_EXCEPTION_FREE_CODE', ...
                                        'SS_OPTION_DISALLOW_CONSTANT_SAMPLE_TIME'};
      % Generate and compile
      aux_sFunMake_casadi( sqpPrepFun.(fieldNames{i}), OptionsCg);
  end

% Add sfunfolder to path
  addpath(genpath(sFunFolder));

fprintf('\nDone\n');

end

