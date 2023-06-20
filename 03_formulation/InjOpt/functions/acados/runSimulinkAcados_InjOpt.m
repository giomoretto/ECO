function uOpt = runSimulinkAcados_InjOpt(ocp,fcn,parModel,parSim,parOpt)

% Prepare parameter struct for simulation
  parECO = getParECO(parOpt,parSim,parModel);
% Get S Functions
  if parOpt.enSFunGen
      getSFun(fcn);
  end
% Get Acados S Functions
  if parOpt.enSFunGen
      getSFunAcados(ocp);
%       getSFunAcados(ocp_single,'single');
  end
  
% Run Simulink simulation
  fprintf('\nRunning Simulink simulation... '); tic
  simOut    = sim('Sim_ECO','SaveOutput','on',...
               'SrcWorkspace','current');
  tSim = toc;fprintf('Done in %.2f seconds \n',tSim);  
  
  % Extract simulation data
  uOpt      = simOut.u.Data(end,:)';
  
end