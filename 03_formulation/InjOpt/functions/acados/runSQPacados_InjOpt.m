function [optVars,status] = runSQPacados_InjOpt(ocp,parModel,parSim,parOpt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    % Define OCP
    ocp = createInitAcadosOCP_InjOpt(ocp,parModel,parSim,parOpt);
    % Call ocp solver
    ocp.solve();
    % Get solution
    xtraj = ocp.get('x');
    status = ocp.get('status');
    switch status
        case 0
            statText = '0 - ''success''';
        case 1
            statText = '1 - ''failure''';
        case 2
            statText = '2 - ''maximum number of SQP iterations reached''';
        case 3
            statText = '3 - ''minimum step size in QP solver reached''';
        case 4
            statText = '4 - ''qp solver failed''';
    end
    ocp.print('stat')
    disp(['Total CPU Time: ',num2str(ocp.get('time_tot')*1000),' ms for ',...
        num2str(ocp.get('sqp_iter')),' SQP Steps. - Status: ',statText]);
    disp(['linearization: ', num2str(ocp.get('time_lin')*1000), ' ms', ...
        ' integrator: ', num2str(ocp.get('time_sim')*1000), ' ms', ...
        ' QP solution: ', num2str(ocp.get('time_qp_sol')*1000), ' ms'])
    % Write acados xTraj into optVars format
    optVarsScaled = xtraj(1:parModel.nInputs,end);
    
    nInj         = parModel.nInputs/2;
    SOENames     = cell(1,nInj); SOENames(1,:) = {'SOE'};
    DOENames     = cell(1,nInj); DOENames(1,:) = {'DOE'};
    uNames       = [SOENames,DOENames];
    optVars      = scaleUnscale(optVarsScaled,uNames,parOpt,'unscale',false);
        
end