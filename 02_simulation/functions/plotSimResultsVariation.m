function [] =  plotSimResultsVariation(fig,simout,parSim)


% Parameters
xlimCA = [-15 35];
% Get Figure
figure(fig);
% Add Data
ax(1) = subplot(2,2,2);grid on;box on; hold on
hPlot(1) = plot(simout.ca,simout.xdot(2,:));
xlabel('[$^\circ \rm CA$]')
ylabel('$Q''_{\rm c}$ [${\rm J}/^\circ \rm CA$]')

ax(2) = subplot(2,2,1);grid on;box on; hold on
hPlot(2) = plot(simout.ca,simout.x(1,:)./1e5);
xlabel('[$^\circ \rm CA$]')
ylabel('$p_{\rm cyl}$ [bar]')

ax(3) = subplot(2,2,3);grid on;box on; hold on
hPlot(3) = plot(simout.ca,simout.y(6,:).*1e6);
title('Ignition Delay')
xlabel('[$^\circ \rm CA$]')
ylabel('$\tau_{\rm ID}$ [$\mu$s]')

ax(4) = subplot(2,2,4);grid on;box on; hold on
colInd = get(ax(4),'ColorOrderIndex');
hPlot(4) = plot(simout.ca,simout.x(2,:),'--'); hold on;
set(ax(4),'ColorOrderIndex',colInd);
hPlot(5) = plot(simout.ca,simout.y(5,:).*42.6e6);
title('Injected vs Released Energy [J]')
xlabel('[$^\circ \rm CA$]')
ylabel('[J]')



linkaxes(ax,'x');
xlim(xlimCA);




end