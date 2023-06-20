function [] = getSFunAcados(ocp,mode)

%% -----------------------------------------------------------------------%
%                    Create ACADOS S-Functions for NLP                    %
%-------------------------------------------------------------------------%

%% Preliminary code
% get build and code folder path, save current location
  buildFolder = what('06_build');
  if exist('mode','var') && strcmp(mode,'single')
      sFunFolder  = fullfile(buildFolder.path,'s_functions','ECO_single');
  else
      sFunFolder  = fullfile(buildFolder.path,'s_functions','ECO');
  end
  codeFolder  = fullfile(sFunFolder,'c_generated_code');
  currentLoc  = pwd;

%% Run S Function Generation
% Print Start Notification
  if exist('mode','var') && strcmp(mode,'single')
      fprintf('\nGenerating Acados S-function (single)... \n\n');
  else
      fprintf('\nGenerating Acados S-function... \n\n');
  end
% Manipulate simulink_opts to [de]activate in- & outputs (preliminary)
  simulink_opts = get_acados_simulink_opts;
  % Inputs
    simulink_opts.inputs.lbx_0          = 1;
    simulink_opts.inputs.ubx_0          = 1;
    simulink_opts.inputs.parameter_traj = 1;
    simulink_opts.inputs.y_ref_0        = 0;
    simulink_opts.inputs.y_ref          = 0;
    simulink_opts.inputs.y_ref_e        = 0;
    simulink_opts.inputs.lbx            = 0;
    simulink_opts.inputs.ubx            = 0;
    simulink_opts.inputs.lbx_e          = 0;
    simulink_opts.inputs.ubx_e          = 0;
    simulink_opts.inputs.lbu            = 0;
    simulink_opts.inputs.ubu            = 0;
    simulink_opts.inputs.lg             = 0;
    simulink_opts.inputs.ug             = 0;
    simulink_opts.inputs.lh             = 1;
    simulink_opts.inputs.uh             = 1;
    simulink_opts.inputs.lh_e           = 1;
    simulink_opts.inputs.uh_e           = 1;
    simulink_opts.inputs.cost_W_0       = 0;
    simulink_opts.inputs.cost_W         = 0;
    simulink_opts.inputs.cost_W_e       = 0;
    simulink_opts.inputs.reset_solver   = 0;
    simulink_opts.inputs.x_init         = 1;
    simulink_opts.inputs.u_init         = 0;
  % outputs
    simulink_opts.outputs.u0            = 0;
    simulink_opts.outputs.utraj         = 0;
    simulink_opts.outputs.xtraj         = 1;
    simulink_opts.outputs.solver_status = 1;
    simulink_opts.outputs.KKT_residual  = 0;
    simulink_opts.outputs.x1            = 0;
    simulink_opts.outputs.CPU_time      = 1;
    simulink_opts.outputs.CPU_time_sim  = 0;
    simulink_opts.outputs.CPU_time_qp   = 0;
    simulink_opts.outputs.CPU_time_lin  = 0;
    simulink_opts.outputs.sqp_iter      = 1;
  % SFun Opts
    simulink_opts.samplingtime = '-1';
    % 't0' (default) - use time step between shooting node 0 and 1
    % '-1' - inherit sampling time from other parts of simulink model
% Go to sFun Folder and generate mex code
  mkdir(sFunFolder);
  cd(sFunFolder);
  ocp.generate_c_code(simulink_opts);
  % ocp.generate_c_code;
% Go to code folder and create s-functions
  cd(codeFolder);
  make_sfun;      % OCP solver
%   make_sfun_sim;  % Integrator
      
% Get back to start folder
  cd(currentLoc);
  
fprintf('\nDone\n');

% Add sfunfolder to path
  addpath(sFunFolder);

% Add codefolder to path
  addpath(codeFolder);

end
