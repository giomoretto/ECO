function plotResults_MATLAB(data)
 
fIdx = 0;

for k = 1:numel(data)
   fIdx      = fIdx + 1;

 % Injector inputs
   SOE1      = data(k).uOpt(1,:)';
   SOE2      = data(k).uOpt(2,:)';
   DOE1      = data(k).uOpt(3,:)';
   DOE2      = data(k).uOpt(4,:)';
   t2ca      = 1e-6*2000/60*2*pi*180/pi;
   SOE4plot  = [SOE1, SOE1+DOE1*t2ca, SOE2, SOE2+DOE2*t2ca];
   u4plot    = repmat([0.7,0,0.7,0],size(SOE4plot,1),1);
   SOE4plot  = [-30*ones(numel(SOE1),1), SOE4plot, 50*ones(numel(SOE1),1)];
   u4plot    = [zeros(numel(SOE1),1), u4plot, zeros(numel(SOE1),1)];

%% Plot Pareto Solutions
   clear s;
   color = {'b','r','g','m'};
   figure(fIdx); clf;
   set(gcf,'Position',[50+(k-1)*50 50+(k-1)*50 560 420]);
   s(1) = subplot('Position',[.09 .65 .69 .27]); hold on; grid on; box on;
   s(1).LineWidth = 1.2;
   s(1).FontSize = 11;
   plot(data(k).y, data(k).eta, 'k--','LineWidth',1.2);
   for kk = 1:numel(data(k).y)
       plot(data(k).y(kk),data(k).eta(kk),'ko','MarkerFaceColor',color{kk},'LineWidth',1.2);
   end
   if k == 2
       set(gca,'xdir','reverse');
   end
   xlim([0.95*min(data(k).y), 1.05*max(data(k).y)]);
   ylim([floor(min(data(k).eta)) ceil(max(data(k).eta))]);
   ylabel('$\eta_{\rm hp}$ [$\%$]');
   xlabel([data(k).name, ' ', data(k).unit]);
   title(['Pareto Solutions for varying ', data(k).name],'FontSize',14);

%% Plot Cylinder Pressure
   s(2) = subplot('Position',[.09 .12 .69 .41]); hold on; grid on; box on;
   s(2).LineWidth = 1.2;
   s(2).FontSize = 11;
   for kk = 1:numel(data(k).y)
       plot(data(k).ca_pCyl, data(k).pCyl(:,kk),color{kk},'LineWidth',1.2);
   end
   axis([data(k).ca_pCyl([1 end]) 30 85]);
   ylabel('$p$ [bar]');
   xlabel('crank angle [$^\circ$CA]');

%% Plot Injector Inputs
   s(3) = subplot('Position',[.78 .72 .2 .2]); hold on; grid on; box on;
   s(3).LineWidth = 1.2;
   s(3).FontSize = 11;
   stairs(SOE4plot(1,:),u4plot(1,:),color{1},'LineWidth',1.2);
   xticklabels('');
   yticklabels('');
   yticks([-1 2]);
   xticks(-20:10:20);
   title('Injector Inputs','FontSize',14)
   
   s(4) = subplot('Position',[.78 .52 .2 .2]); hold on; grid on; box on;
   s(4).LineWidth = 1.2;
   s(4).FontSize = 11;
   stairs(SOE4plot(2,:),u4plot(2,:),color{2},'LineWidth',1.2);
   xticklabels('');
   yticklabels('');
   yticks([-1 2]);
   xticks(-20:10:20);

   s(5) = subplot('Position',[.78 .32 .2 .2]); hold on; grid on; box on;
   s(5).LineWidth = 1.2;
   s(5).FontSize = 11;
   stairs(SOE4plot(3,:),u4plot(3,:),color{3},'LineWidth',1.2);
   xticklabels('');
   yticklabels('');
   yticks([-1 2]);
   xticks(-20:10:20);

   s(6) = subplot('Position',[.78 .12 .2 .2]); hold on; grid on; box on;
   s(6).LineWidth = 1.2;
   s(6).FontSize = 11;
   stairs(SOE4plot(4,:),u4plot(4,:),color{4},'LineWidth',1.2);
   stairs([-50 50],[0 0],'k','LineWidth',1.2);
   yticklabels('');
   xlabel('crank angle [$^\circ$CA]');
   yticks([-1 2]);
   xticks(-20:10:20);

   linkaxes(s(3:6))
   axis([-17 20 0 1]);
   
end