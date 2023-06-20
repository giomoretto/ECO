function simout = CompleteSimulation(ca,x0,u,parSim,parModel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% global test

 % Initialize Vectors
   x      = zeros(parModel.nStates,numel(ca));
   y      = zeros(parModel.nOutputs,numel(ca));
   xdot   = zeros(parModel.nStates,numel(ca));
   x(:,1) = x0;
 % Run Simulation
   for kk = 1:numel(ca)
       xNew     = x(:,kk);
       caAct    = ca(kk);
       % Get current deltaPhi
       if kk ~= numel(ca)
           deltaPhi = (ca(kk+1)-ca(kk))/parSim.Opts.nInt;
       else
           deltaPhi = (ca(kk)-ca(kk-1))/parSim.Opts.nInt;
       end
       % Do integration
       for jj = 1:parSim.Opts.nInt
           switch parSim.Opts.Int
               case 'EulerFW'
                   if jj == 1
                       [xNew,xdot(:,kk),y(:,kk)]  = ...
                           EulerFWIntegrationStep(xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   else
                       [xNew,~,~]  = EulerFWIntegrationStep(...
                           xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   end
% 
% test.RK.uNOx(test.j) = test.uNOx(end-3);
%                    test.RK.NOscale(test.j) = test.NOscale(end);
%                    test.RK.ca(test.j) = caAct;
%                    test.RK.Tbz(test.j) = test.thetaBZ(end);
%                    test.RK.Vbz(test.j) = test.Vbz(end);
%                    test.RK.dVbz(test.j) = test.dVbz(end);
%                    test.RK.dNOchem(test.j) = test.dNOchem(end);
%                    test.RK.dNOdV(test.j) = test.dNOdV(end);
%                    test.RK.dNOdphi(test.j) = test.dNOdphi(end);
%                    test.RK.zetaComb(test.j) = test.zetaComb(end);
%                    test.RK.xBZ(test.j) = test.xBZ(end);
%                    test.RK.xBZsat(test.j) = test.xBZsat(end);
%                    test.RK.xBG(test.j) = test.xBG(end);
%                    test.RK.mf2c(test.j) = test.mf2c(end);
% 
% 
%                    test.j = test.j + 1;

               case 'RK4'
                   if jj == 1                  

                       [xNew,xdot(:,kk),y(:,kk)]  = RK4IntegrationStep(...
                           xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   else
                       [xNew,~,~]  = RK4IntegrationStep(...
                           xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   end
                
%                    test.RK.uNOx(test.j) = test.uNOx(end-3);
%                    test.RK.NOscale(test.j) = test.NOscale(end-3);
%                    test.RK.ca(test.j) = caAct;
%                    test.RK.Tbz(test.j) = test.thetaBZ(end-3);
%                    test.RK.Vbz(test.j) = test.Vbz(end-3);
%                    test.RK.dVbz(test.j) = test.dVbz(end-3);
%                    test.RK.dNOchem(test.j) = test.dNOchem(end-3);
%                    test.RK.dNOdV(test.j) = test.dNOdV(end-3);
%                    test.RK.dNOdphi(test.j) = test.dNOdphi(end-3);
%                    test.RK.zetaComb(test.j) = test.zetaComb(end-3);
%                    test.RK.xBZ(test.j) = test.xBZ(end-3);
%                    test.RK.xBZsat(test.j) = test.xBZsat(end-3);
%                    test.RK.xBG(test.j) = test.xBG(end-3);
%                    test.RK.mf2c(test.j) = test.mf2c(end-3);
%                    test.RK.nTot(test.j) = test.nTot(end-3);
%                    test.RK.R1conc(test.j) = test.R1conc(end-3);
%                    test.RK.R2conc(test.j) = test.R2conc(end-3);
%                    test.RK.R3conc(test.j) = test.R3conc(end-3);
%                    test.RK.Qconc(test.j) = test.Qconc(end-3);
%                    test.RK.NOeq(test.j) = test.NOeq(end-3);
%                    test.RK.X_NO(test.j) = test.X_NO(end-3);
%                    test.RK.k1f(test.j) = test.k1f(end-3);
%                    test.RK.X_N2(test.j) = test.X_N2(end-3);
%                    test.RK.X_O(test.j) = test.X_O(end-3);
% 
% 
%                    test.j = test.j + 1;


               case 1 % Simulink implementation for RK4
                   if jj == 1
                       [xNew,xdot(:,kk),y(:,kk)]  = RK4IntegrationStep(...
                           xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   else
                       [xNew,~,~]  = RK4IntegrationStep(...
                           xNew,u(:,kk),caAct,...
                           parModel,parSim.OP,deltaPhi);
                   end
           end
           caAct = caAct + deltaPhi;
       end
       if kk ~= numel(ca)
           x(:,kk+1) = xNew;
       end
   end
 % Write Results to struct
   simout.ca   = ca;
   simout.x    = x;
   simout.u    = u;
   simout.y    = y;
   simout.xdot = xdot;
end
%
%
function [xNew,xdot,y] = EulerFWIntegrationStep(xAct,uAct,caAct,parModel,parOP,deltaPhi)
    [xdot,y] = CompleteModel(xAct,uAct,caAct,...
                               parModel,parOP,'Num');
    xNew = xAct+xdot*deltaPhi;
end
%
%
function [xNew,k1,y] = RK4IntegrationStep(xAct,uAct,caAct,parModel,parOP,deltaPhi)
     % Define k1
       [k1,y] = CompleteModel(xAct,uAct,caAct,parModel,parOP,'Num');
     % Calculate k2
       x_k1 = xAct+k1*deltaPhi/2;
       [k2,~] = CompleteModel(x_k1,uAct,caAct+deltaPhi/2,parModel,parOP,'Num');
     % Calculate k3
       x_k2 = xAct+k2*deltaPhi/2;
       [k3,~] = CompleteModel(x_k2,uAct,caAct+deltaPhi/2,parModel,parOP,'Num');
     % Calculate k4
       x_k3 = xAct+k3*deltaPhi;
       [k4,~] = CompleteModel(x_k3,uAct,caAct+deltaPhi,parModel,parOP,'Num');
     % Calculate New x
       xNew = xAct + deltaPhi/6 * (k1 + 2*k2 + 2*k3 + k4);
    
end
