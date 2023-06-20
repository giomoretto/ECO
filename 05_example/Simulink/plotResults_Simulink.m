function plotResults_Simulink(ref, simOut, parECO)

%% Extract and Prepare Simulation Data 
 % Reference values
   % Remove initial sequence necessary for converged IC
   ref.imep  = [ref.imep(3:end),  ref.imep(end) ] * 1e-5;
   ref.dpMax = [ref.dpMax(3:end), ref.dpMax(end)] * 1e-5;
   ref.Tevo  = [ref.Tevo(3:end),  ref.Tevo(end) ] - 273;
   ref.cNOx  = [ref.cNOx(3:end),  ref.cNOx(end) ];
   % Reorder for visualization
   ref.time  = [0, repelem(ref.time(3:end),2), 0];
   ref.imep  = [repelem(ref.imep(1:end-1),2), -15, -15];
   ref.dpMax = [repelem(ref.dpMax(1:end-1),2), 20, 20];
   ref.Tevo  = [repelem(ref.Tevo(1:end-1),2), 0, 0];
   ref.cNOx  = [repelem(ref.cNOx(1:end-1),2), 1e4, 1e4];
 % Simulation outputs
   idx_y	 = find(simOut.y.Time >= ref.dT,1,'first');
   y.time    = simOut.y.Time(1:end-idx_y);
   y.imep    = simOut.y.Data(idx_y+1:end,1) * 1e-5;
   y.dpMax   = simOut.y.Data(idx_y+1:end,4) * 1e-5;
   y.Tevo    = simOut.y.Data(idx_y+1:end,3) - 273;
   y.cNOx    = simOut.y.Data(idx_y+1:end,2);
 % Injector inputs
   u.t2ca    = 1e-6*2000/60*2*pi*180/pi;
   u.SOE1    = simOut.u.Data(idx_y+1:end,1);
   u.SOE2    = simOut.u.Data(idx_y+1:end,2);
   u.DOE1    = simOut.u.Data(idx_y+1:end,3);
   u.DOE2    = simOut.u.Data(idx_y+1:end,4);
   u.SOE     = [u.SOE1, u.SOE1+u.DOE1*u.t2ca, u.SOE2, u.SOE2+u.DOE2*u.t2ca];
   u.y       = repmat([0.7,0,0.7,0],size(u.SOE,1),1);
   u.y(u.DOE1 < 101,1) = 0;
   u.SOE     = [-30*ones(numel(u.SOE1),1), u.SOE, 50*ones(numel(u.SOE1),1)];
   u.y       = [zeros(numel(u.SOE1),1), u.y, zeros(numel(u.SOE1),1)];
 % Simulation for cylinder pressure
   y.caSim   = -172:0.5:30;
   y.x0Sim   = [parECO.parSim.OP.pInt,0,0,parECO.parSim.OP.thetaIvc,0]';
   
%% Prepare Plot Appearance
   pos.f		= [284 257 869 475];
   pos.scale	= pos.f(3)/pos.f(4);
   pos.x0		= [0.08 0.7];
   pos.y0		= ones(1,2)*pos.scale*pos.x0(1);
   pos.x		= [0.55 0.27];
   pos.y		= [0.2 0.55 0.25];

%% Create GIF
   delta = 50;
   figure('Position',pos.f,'Visible','off');
   
 % Create figures
   F = struct('cdata',cell(1, numel(y.time)/delta),'colormap',cell(1, numel(y.time)/delta));
   fprintf('Prepare Visualization.. ');
   for idx = [1:delta:numel(y.time) numel(y.time)]
       clf;
     % Plot imep
       y.imepTemp = nan(size(y.time));
       y.imepTemp(1:idx) = y.imep(1:idx);
       axes('Position',[pos.x0(1) pos.y0(1)+3*pos.y(1) pos.x(1) pos.y(1)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       patch(ref.time, ref.imep,[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.3,'HandleVisibility','off')
       plot(y.time,y.imepTemp,'k','LineWidth',1.2);
       axis([ref.time([1 end-1]) min(y.imep)-0.4 max(y.imep)+0.4])
       xticklabels('');
       ylabel('$w_{\rm hp}$ [bar]');
     % Plot dpMax
       y.dpMaxTemp = nan(size(y.time));
       y.dpMaxTemp(1:idx) = y.dpMax(1:idx);
       axes('Position',[pos.x0(1) pos.y0(1)+2*pos.y(1) pos.x(1) pos.y(1)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       patch(ref.time, ref.dpMax,[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.3,'HandleVisibility','off')
       plot(y.time,y.dpMaxTemp,'k','LineWidth',1.2);
       axis([ref.time([1 end-1]) min(y.dpMax)-0.4 max(y.dpMax)+0.4])
       xticklabels('');
       ylabel('${\rm max}(p'')$ [$\frac{\rm bar}{\rm ^\circ CA}$]')
     % Plot Tevo
       y.TevoTemp = nan(size(y.time));
       y.TevoTemp(1:idx) = y.Tevo(1:idx);
       axes('Position',[pos.x0(1) pos.y0(1)+1*pos.y(1) pos.x(1) pos.y(1)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       patch(ref.time, ref.Tevo,[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.3,'HandleVisibility','off')
       plot(y.time,y.TevoTemp,'k','LineWidth',1.2);
       axis([ref.time([1 end-1]) min(y.Tevo)-30 max(y.Tevo)+30])
       xticklabels('');
       ylabel('$T_{\rm evo}$ [$^\circ$C]')
     % Plot cNOx
       y.cNOxTemp = nan(size(y.time));
       y.cNOxTemp(1:idx) = y.cNOx(1:idx);
       axes('Position',[pos.x0(1) pos.y0(1) pos.x(1) pos.y(1)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       patch(ref.time, ref.cNOx,[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeAlpha',0.3,'HandleVisibility','off')
       plot(y.time,y.cNOxTemp,'k','LineWidth',1.2);
       axis([ref.time([1 end-1]) min(y.cNOx)-100 max(y.cNOx)+100])
       xlabel('time [s]');
       ylabel('$X_{\rm NO}$ [ppm]')
     % Plot pCyl
	   y.uSim = repmat([u.SOE1(idx);u.SOE2(idx);u.DOE1(idx);u.DOE2(idx)],1,numel(y.caSim));
	   simOut = CompleteSimulation(y.caSim,y.x0Sim,y.uSim,parECO.parSim,parECO.parModel);
       axes('Position',[pos.x0(2) pos.y0(2)+pos.y(3) pos.x(2) pos.y(2)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       plot(-20:0.5:30,1e-5*simOut.x(1,simOut.ca >= -20 & simOut.ca <= 30),'k','LineWidth',1.2);
       axis([-20 30 20 100]);
       xticklabels('');
       ylabel('$p$ [bar]')
     % Plot injector inputs 
       axes('Position',[pos.x0(2) pos.y0(1) pos.x(2) pos.y(3)]); hold on; grid on; box on;
       set(gca,'LineWidth',1.2,'FontSize',11);
       stairs(u.SOE(idx,:),u.y(idx,:),'k','LineWidth',1.5)
       axis([-20 30 0 1]);
       yticklabels('');
       xlabel('crank angle [$^\circ$CA]')
     % Save current frame
       set(gcf,'Color','w')
       F(floor(idx/delta+1)) = getframe(gcf);
   end
   close(gcf);
   fprintf('Done!\n\n');
   
 % Play Movie
   fprintf('Configure and Run Visualization.. ');
   movie = Movie(F);
   movie.PlayMovie();
   fprintf('\n\nVisualization closed.\n\n');
 
end