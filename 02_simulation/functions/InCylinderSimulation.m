function simout = InCylinderSimulation(ca,x0,dQComb,parSim,parModel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

 % Initialize Vectors
   u      = dQComb;
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
    [xdot,y] = InCylinderModel(xAct,uAct,caAct,...
                               parModel,parOP,'Num');
    xNew = xAct+xdot*deltaPhi;
end
%
%
function [xNew,k1,y] = RK4IntegrationStep(xAct,uAct,caAct,parModel,parOP,deltaPhi)
     % Define k1
       [k1,y] = InCylinderModel(xAct,uAct,caAct,parModel,parOP,'Num');
     % Calculate k2
       x_k1 = xAct+k1*deltaPhi/2;
       [k2,~] = InCylinderModel(x_k1,uAct,caAct+ ...
           deltaPhi/2,parModel,parOP,'Num');
     % Calculate k3
       x_k2 = xAct+k2*deltaPhi/2;
       [k3,~] = InCylinderModel(x_k2,uAct,caAct+ ...
           deltaPhi/2,parModel,parOP,'Num');
     % Calculate k4
       x_k3 = xAct+k3*deltaPhi;
       [k4,~] = InCylinderModel(x_k3,uAct,caAct+ ...
           deltaPhi,parModel,parOP,'Num');
     % Calculate New x
       xNew = xAct + deltaPhi/6 * (k1 + 2*k2 + 2*k3 + k4);
end
