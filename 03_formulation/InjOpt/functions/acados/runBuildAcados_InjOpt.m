function [] = runBuildAcados_InjOpt(ocp,fcn,parModel,parSim,parOpt)

% Get build folder path
  buildFolder	= what('06_build');
  sFunFolder	= fullfile(buildFolder.path,'s_functions');

% Prepare parameter struct for simulation
  parECO = getParECO(parOpt,parSim,parModel);
  save(fullfile(buildFolder.path,'parECO.mat'),'parECO');
  
% Get S Functions
  if parOpt.enSFunGen
      getSFun(fcn);
  end
% Get Acados S Functions
  if parOpt.enSFunGen
      getSFunAcados(ocp);
  end
  
% Save parECO.mat and s-functions
  fileName	= '';
  if ~isempty(fileName)
	  copyfile(sFunFolder,fullfile(buildFolder.path,['s_functions_',fileName]));
	  copyfile(fullfile(buildFolder.path,'parECO.mat'), ...
		  fullfile(buildFolder.path,['parECO_',fileName,'.mat']));
  end
  
end