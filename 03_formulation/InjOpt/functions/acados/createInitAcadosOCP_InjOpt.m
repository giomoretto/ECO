function [ocp] = createInitAcadosOCP_InjOpt(ocp,parModel,parSim,parOpt)
% Calculate x0 and initial trajectory with simulation
  % Prepare Simulation
    caPreOpt     = (parSim.OP.ca_ivc:parSim.Opts.deltaPhi:...
                        parOpt.OptimizationRange(1));
    u0PreSim     = repmat(parOpt.u0,1,numel(caPreOpt));
    x0PreOpt     = [parSim.OP.pInt,0,0]';			% [p QComb IMEP]
    if parOpt.enNOx
        x0PreOpt = [x0PreOpt;parSim.OP.thetaIvc;0]; % [p QComb IMEP thetaUZ NOx]
    end
  % Run simulation
    simPreOpt    = CompleteSimulation(caPreOpt,x0PreOpt,u0PreSim,parSim,parModel);
  % Get x0 at Start of Opt Range
    x0Opt        = simPreOpt.x(:,end);
  % Calculate states for init trajectory with simulation
    u0Sim        = repmat(parOpt.u0,1,numel(parOpt.ca));
    simOpt0      = CompleteSimulation(parOpt.ca,x0Opt,u0Sim,parSim,parModel);
% Define and scale needed parameters
  % States
    if parOpt.enNOx
        xNames   = {'pCyl','QComb','IMEP','Theta','NO'};
    else
        xNames   = {'pCyl','QComb','IMEP'};
    end
    [x0Opt,xScale,xOffs] = scaleUnscale(x0Opt,xNames,parOpt,'scale',false);
    xScaleVec    = repmat(xScale,1,size(simOpt0.x,2));
    xOffsVec     = repmat(xOffs,1,size(simOpt0.x,2));
    x_traj_init  = (simOpt0.x-xOffsVec)./xScaleVec;
  % Inputs
    nInj         = parModel.nInputs/2;
    SOENames     = cell(1,nInj); SOENames(1,:) = {'SOE'};
    DOENames     = cell(1,nInj); DOENames(1,:) = {'DOE'};
    uNames       = [SOENames,DOENames];
    [~,uScale,uOffs] = scaleUnscale(parOpt.u0,uNames,parOpt,'scale',false);
    uScaleVec    = repmat(uScale,1,size(u0Sim,2));
    uOffsVec     = repmat(uOffs,1,size(u0Sim,2));
    u_traj_init  = (u0Sim-uOffsVec)./uScaleVec;
  % Crank-angle
    [~,caScale,caOffs] = scaleUnscale(parOpt.OptimizationRange(1),{'SOE'},parOpt,'scale',false);
    caScaleVec   = repmat(caScale,1,numel(parOpt.ca));
    caOffsVec    = repmat(caOffs,1,numel(parOpt.ca));
    ca_traj_init = (parOpt.ca-caOffsVec)./caScaleVec;
% Update constraints
  % Prepare constraints
    infRep		 = 1e2;
    N			 = ocp.opts_struct.param_scheme_N;
    % Actual states
    lbx			 = scaleUnscale(parOpt.xMin,xNames,parOpt,'scale',false);
    ubx			 = scaleUnscale(parOpt.xMax,xNames,parOpt,'scale',false);
    % Inputs
    lbu			 = scaleUnscale(parOpt.uMin,uNames,parOpt,'scale',false);
    ubu			 = scaleUnscale(parOpt.uMax,uNames,parOpt,'scale',false);
    % Crank-angle
    calb		 = scaleUnscale(parOpt.ca(1),{'SOE'},parOpt,'scale',false);
    caub		 = scaleUnscale(parOpt.ca(1),{'SOE'},parOpt,'scale',false);
    % References
    imep		 = scaleUnscale(parOpt.Reference.IMEP,{'IMEP'},parOpt,'scale',false);
    Tevo		 = scaleUnscale(parOpt.Reference.Tmin,{'Theta'},parOpt,'scale',false);
    pMax		 = scaleUnscale(parOpt.Reference.pMax,{'pCyl'},parOpt,'scale',false);
    dpMax		 = scaleUnscale(parOpt.Reference.dpMax,{'pCyl'},parOpt,'scale',true);
    cNOx		 = scaleUnscale(parOpt.Reference.cNOx,{'NOppm'},parOpt,'scale',false);
    PhiMax		 = parOpt.Reference.PhiMax;
    % Latest center uf combustion
    idxCoC		 = find(parOpt.ca < parOpt.Reference.CoCmax,1,'last');
  % Write to model object: constraints
    % Initial constraints (consisting of actual states and inputs)
    ocp.set('constr_lbx', [lbu; x0Opt; calb], 0);
    ocp.set('constr_ubx', [ubu; x0Opt; caub], 0);
    % State constraints (consisting of actual states and inputs)
    % for ii = 1:N-1
    %     ocp.set('constr_lbx', lbx, ii);
    %     ocp.set('constr_ubx', ubx, ii);
    % end
    % Reference constraints
    % pMax, dpMax and CoC
    for ii = 1:N
		ocp.set('constr_uh', [pMax; dpMax; infRep], ii-1);
        ocp.set('constr_lh', [0; -infRep; -infRep], ii-1);
        if ii > idxCoC
            ocp.set('constr_lh', [0; -infRep; 0], ii-1);
        end
    end
    % imep, Tevo, (NOx), Phi and injection distancing
    if parOpt.enNOx
        ocp.set('constr_lh', [imep; Tevo; 0; 0; zeros(nInj-1,1)], N);
        ocp.set('constr_uh', [infRep; infRep; cNOx; PhiMax; infRep*ones(nInj-1,1)], N);
    else
        ocp.set('constr_lh', [imep; Tevo; 0; zeros(nInj-1,1)], N);
    end
% Set trajectory initialization for very first ocp
  ux_traj_init  = [u_traj_init; x_traj_init; ca_traj_init];
  if isnan(ocp.get_cost)
      ocp.set('init_x', ux_traj_init);
  end
% Set parameters (OP and dt_Inj)
  p           = [parSim.OP.engSpd;			  ...
                 parSim.OP.pIM;				  ...
                 parSim.OP.thetaIM;			  ...
                 parSim.OP.xiBgIM;			  ...
                 parSim.OP.pRail;			  ...
                 parModel.wall.scaleFiredWHL; ...
				 parOpt.dt_Inj];
  for ii = 1:N+1
      ocp.set('p',p,ii-1);
  end
end