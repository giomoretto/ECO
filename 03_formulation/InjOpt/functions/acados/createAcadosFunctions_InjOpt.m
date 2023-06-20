function [ocp,fcn] = createAcadosFunctions_InjOpt(parOpt,parModel)

% Initialization
  import casadi.*
  nStates          = parModel.nStates;
  nInputs          = parModel.nInputs;
  nInj             = nInputs/2;
  caNum            = parOpt.ca;
  N                = numel(caNum)-1; % last node is simulated and constrained with terminal constraint
  T                = sum(diff(caNum));  

% Define acados settings
  if parOpt.SQP.nSQPmax == 1
      nlp_solver   = 'sqp_rti';
  else
      nlp_solver   = 'sqp';
  end
  qp_solver        = 'partial_condensing_hpipm';
  qp_solver_cond_N = 5;         % For partial condensing
  model_name       = 'ECO';
  sim_method       = 'erk';     % erk, irk, irk_gnsf
  sim_method_num_stages = 4;    % Runge-Kutta int. stages: (1) RK1, (4) RK4
  cost_type        = 'ext_cost';

%% Create Acados Model Object
% Initialize acados model object
  ocp_model     = acados_ocp_model();
  ocp_model.set('name', model_name);
  ocp_model.set('T', T);

% Create symbolic model expressions with casADi
  % States (consisting of actual states and inputs, scaled)
    uxSym       = SX.sym('x',nInputs+nStates+1,1);
    uxdotSym    = SX.sym('xdot',nInputs+nStates+1,1);
  % Parameters (Operating Point and Injection Distancing)
    p           = SX.sym('p',5+2,1);
  % Extract inputs and actual states
    u           = uxSym(1:nInputs);
    x           = uxSym(nInputs+1:nInputs+nStates);
    ca          = uxSym(nInputs+nStates+1);
  % Extract OP and recalculate intake conditions
    OP.Ne       = p(1);
    OP.p_im     = p(2);
    OP.T_im     = p(3);
    OP.x_bg     = p(4);
    OP.p_rail   = p(5);
    parOP       = parOPDef(parModel,OP);
  % Extract scaling for WHL
    parModel.wall.scaleFiredWHL = p(5+1);
  % Extract injection distancing
    dtInj       = p(5+2);
  % Unscale states and inputs before sending them to model function
    if parOpt.enNOx
        xNames  = {'pCyl','QComb','IMEP','Theta','NO'};
        yNames  = {'Theta','kappa','dQWall','specR','mFuelPrep','tauIgn','xBG','Theta','NOppm'};
    else
        xNames  = {'pCyl','QComb','IMEP'};
        yNames  = {'Theta','kappa','dQWall','specR','mFuelPrep','tauIgn'};
    end
    SOENames    = cell(1,nInj); SOENames(1,:) = {'SOE'};
    DOENames    = cell(1,nInj); DOENames(1,:) = {'DOE'};
    xUnscaled   = scaleUnscale(x,xNames,parOpt,'unscale',false);
    SOEUnscaled = scaleUnscale(u(1:nInj),SOENames,parOpt,...
                    'unscale',false);
    DOEUnscaled = scaleUnscale(u(nInj+1:end),DOENames,parOpt,...
                    'unscale',false);
    uUnscaled   = [SOEUnscaled;DOEUnscaled];
    caUnscaled  = scaleUnscale(ca,{'SOE'},parOpt,'unscale',false);
  % System dynamics
    [xdotUnscaled,yUnscaled]  = CompleteModel(...
                    xUnscaled,uUnscaled,caUnscaled,parModel,parOP,'Sym');
    xdot        = scaleUnscale(xdotUnscaled,xNames,parOpt,'scale',true);
    y           = scaleUnscale(yUnscaled,yNames,parOpt,'scale',false);
  % Get explicit and implicit derivative expressions for ux vector
    uxdot       = [zeros(nInputs,1); xdot; scaleUnscale(1,{'SOE'},parOpt,'scale',true)];
    uxdot_impl  = uxdot - uxdotSym;
  % Write to model object: symbolics
    ocp_model.set('sym_x', uxSym);
    ocp_model.set('sym_xdot', uxdotSym);
    ocp_model.set('sym_p',p);
  % Write to model object: dynamics
    if (strcmp(sim_method, 'erk'))
        ocp_model.set('dyn_type', 'explicit');
        ocp_model.set('dyn_expr_f', uxdot);
    else % irk irk_gnsf
        ocp_model.set('dyn_type', 'implicit');
        ocp_model.set('dyn_expr_f', uxdot_impl);
    end
    
% Create symbolic expression for constraints    
  % Prepare constraints
    infRep      = 1e2;
    % Actual states
    lbx         = -infRep * ones(size(parOpt.xMin));
    ubx         =  infRep * ones(size(parOpt.xMax));
    Jbx         = [zeros(numel(x),numel(u)), eye(numel(x))];
    % Inputs
    lbu         = -infRep * ones(size(parOpt.uMin));
    ubu         =  infRep * ones(size(parOpt.uMax));
    % States (consisting of actual states and inputs), for initial condition
    lbux        = [lbu; lbx; -infRep];
    ubux        = [ubu; ubx; infRep];
    Jbux        = eye(numel(uxSym));
    % Injection distancing
    bInj        = [];
    [~,mtot,DOI_mus] = AlgebraicInjectorModel(0,parOP.pRail,SOEUnscaled,DOEUnscaled);
    if nInj > 1
        % Get min. distance between Injections scaled
        dCA_Inj = (DOI_mus+dtInj) * parOP.engSpd*360*1e-6;    % degCA
        dCA_Inj = scaleUnscale(dCA_Inj,SOENames,parOpt,'scale',true); % scaled
        % Iterate through Injections and add constraint
        for i = 2:nInj
            bInj      = [bInj;u(i)-u(i-1)-dCA_Inj(i-1)]; %#ok<AGROW>
        end
    end
    % Latest center of combustion
    Qtot        = mtot*parModel.thermo.fuel.lowHeatVal;
    h_coc       = x{2}-Qtot/2/parOpt.Scale(2);
    % calculate FA-ratio
    Phi = mtot*parModel.thermo.fuel.airFuelRatioSt*parModel.thermo.xi.air.O2...
        / parOP.mCylTot/parOP.xiO2;
    
  % Write to model object: constraints
    % Initial constraints (consisting of actual states and inputs)
    ocp_model.set('constr_lbx_0', lbux);
    ocp_model.set('constr_ubx_0', ubux);
    ocp_model.set('constr_Jbx_0', Jbux);
    % State constraints (only actual states)
    % ocp_model.set('constr_lbx', lbx);
    % ocp_model.set('constr_ubx', ubx);
    % ocp_model.set('constr_Jbx', Jbx);
    % Reference constraints
    % pMax, dpMax and CoC
    ocp_model.set('constr_type', 'bgh');
    ocp_model.set('constr_expr_h',[x{1}; xdot{1}; h_coc]);
    ocp_model.set('constr_lh', -infRep*ones(3,1));
    ocp_model.set('constr_uh', infRep*ones(3,1));
    % imep, Tevo, (NOx) and injection distancing
    if parOpt.enNOx
        ocp_model.set('constr_type_e', 'bgh');
        ocp_model.set('constr_expr_h_e',[x{3}; y{1}; y{9}; Phi; bInj]);
        ocp_model.set('constr_lh_e', -infRep*ones(4 + nInj-1,1));
        ocp_model.set('constr_uh_e', infRep*ones(4 + nInj-1,1));
    else
        ocp_model.set('constr_type_e', 'bgh');
        ocp_model.set('constr_expr_h_e',[x{3}; y{1}; Phi; bInj]);
        ocp_model.set('constr_lh_e', -infRep*ones(3 + nInj-1,1));
        ocp_model.set('constr_uh_e', infRep*ones(3 + nInj-1,1));
    end
    % Slack on NOx
    if parOpt.enNOx
        ocp_model.set('constr_Jsh_e',[0; 0; 1; 0; zeros(nInj-1,1)]);
        ocp_model.set('cost_Zl_e',1e3);
        ocp_model.set('cost_Zu_e',10);
        ocp_model.set('cost_zl_e',1e3);
        ocp_model.set('cost_zu_e',10);
    end
    % Slack on dpMax
    ocp_model.set('constr_Jsh',[0; 1; 0]);
    ocp_model.set('cost_Zl',1e3);
    ocp_model.set('cost_Zu',20);
    ocp_model.set('cost_zl',1e3);
    ocp_model.set('cost_zu',1);
% % ... see ocp_model.model_struct to see what other fields can be set

% Create Initial Variables for Cost
  % Define cost function
    cost_e      = scaleUnscale(Qtot,{'QComb'},parOpt,'scale',false);
  % Write to model object: cost
    ocp_model.set('cost_type', cost_type);
    ocp_model.set('cost_type_e', cost_type);
    ocp_model.set('cost_expr_ext_cost', 0);        % Stage cost
    ocp_model.set('cost_expr_ext_cost_e',cost_e);  % Terminal cost
  
%% Create Acados Option Object
% Std Settings
  ocp_opts = acados_ocp_opts(); 

  ocp_opts.set('param_scheme_N', N);
  if any(diff(caNum) ~= parOpt.deltaPhi)
      ocp_opts.set('time_steps',diff(caNum));
  end
  ocp_opts.set('nlp_solver', nlp_solver);
  ocp_opts.set('sim_method', sim_method);
  ocp_opts.set('sim_method_num_stages',sim_method_num_stages);
  ocp_opts.set('sim_method_num_steps', parOpt.nInt);
  ocp_opts.set('qp_solver', qp_solver);
  ocp_opts.set('qp_solver_cond_N', qp_solver_cond_N);
  ocp_opts.set('nlp_solver_max_iter',parOpt.SQP.nSQPmax);

% SQP Step size & Globlalization
  ocp_opts.set('levenberg_marquardt', 1e-2);
  ocp_opts.set('nlp_solver_step_length',parOpt.SQP.StepSize );

% Tolerances
  tolOpts = {'nlp_solver_tol_stat','nlp_solver_tol_eq',...
             'nlp_solver_tol_ineq','nlp_solver_tol_comp'};
  for i =1:numel(tolOpts)
      switch i
          case 1
              ocp_opts.set(tolOpts{i},1e-4);
          otherwise
              ocp_opts.set(tolOpts{i},1e-6);
      end
  end
        
% Set output directory
  buildFolder = what('06_build');
  ocp_opts.set('output_dir',fullfile(buildFolder.path,'acadosBuild','ECO'));
  % ... see ocp_opts.opts_struct to see what other fields can be set
    
%% Create Acados OCP Object
% Create OCP
  ocp			= acados_ocp(ocp_model, ocp_opts);
  ocp.opts_struct.Tsim = parOpt.deltaPhi;
  
%% Create Casadi Functions (for simulink simulation)
% Create casadi functions for model dynamics
  fuxdot		= Function('fuxdot',{uxSym,p(1:6)},{uxdot});
  fuxdot		= fuxdot.expand();
% Create casadi functions for model dynamics
  fy			= Function('fy',{uxSym,p(1:6)},{y});

% Create casadi functions for model simulation
  fsim			= createAcadosSimulation(fuxdot,fy,parOpt,parModel);
  fcn.fsim		= fsim.expand();

end