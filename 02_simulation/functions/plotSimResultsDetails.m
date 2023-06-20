function [] = plotSimResultsDetails(fig3,simOpt,parSim,parOpt)

set(fig3,'Position',[745.0000  146.1429  557.7143  581.7143])

nInputs = size(simOpt.u,1)/2;

% plot active injections
t2ca = parSim.OP.engSpd*2*pi*180/pi;
uSOE = simOpt.u(1:nInputs,1);
uDOE = simOpt.u(nInputs+1:end,1)*1e-6;
boolU = uDOE(:,1) < (parOpt.uMin(end)+1)*1e-6;

uSOE(boolU,:) = [];
uDOE(boolU,:) = [];

for i = 1:numel(uSOE(:,1))
    uSOE4plot(2*i-1) = uSOE(i,1);
    uSOE4plot(2*i)   = uSOE(i,1)+uDOE(i,1)*t2ca;
    u4plot(2*i-1) = 0.3;
    u4plot(2*i) = 0;
end

uSOE4plot = [-30, uSOE4plot];
u4plot = [0, u4plot];

% simulation resuts
pMax  = max(simOpt.x(1,:)*1e-5);
dpMax = max(simOpt.xdot(1,:)*1e-5);
[~,idxCA50] = min((simOpt.x(2,:)-simOpt.x(2,end)/2).^2);
CA50 = simOpt.ca(idxCA50);
imep = max(simOpt.x(3,end))*1e-5;
if parOpt.enNOx
    Tbz  = max(simOpt.y(8,:));
    NOx  = simOpt.y(9,end);
    Phi  = simOpt.y(10,end);
else
    Phi = simOpt.y(7,end);
end

% plot results 
mFig = 3;
nFig = 2;

s(1) = subplot(mFig,nFig,1); hold on; grid on; box on
plot(simOpt.ca,  simOpt.x(1,:)*1e-5,'k')
xlabel('[$^\circ \rm CA$]')
ylabel('$p_{\rm cyl}$ [bar]')
yyaxis right
plot(simOpt.ca,  simOpt.xdot(1,:)*1e-5)
ylabel('$p''_{\rm cyl}$ [${\rm bar}/^\circ \rm CA$]')

s(2) = subplot(mFig,nFig,3); hold on; grid on; box on
plot(simOpt.ca,  simOpt.xdot(2,:),'k')
ylabel('$Q''_{\rm c}$ [${\rm J}/^\circ \rm CA$]')
xlabel('[$^\circ \rm CA$]')


yyaxis right
stairs(uSOE4plot,u4plot,'k--')
ylim([0 1])


hAxes = get(gca);
set(hAxes.YAxis(2),'visible','off')


s(3) = subplot(mFig,nFig,5); hold on; grid on; box on
plot(simOpt.ca, simOpt.y(1,:),'k')
xlabel('[$^\circ \rm CA$]')
ylabel('T [K]')

if parOpt.enNOx
    plot(simOpt.ca, simOpt.x(4,:),'b')
    plot(simOpt.ca, simOpt.y(8,:),'r')

    yyaxis right
    plot(simOpt.ca, simOpt.y(9,:),'g')
    ylabel('$X_{\rm NO}$ [ppm]')

end


linkaxes(s,'x')
xlim([-15 35])


strIdx = 1;
strPlt{strIdx} = '\textbf{Model Parameters, Operating Point}';

for i = 2:5
    strIdx = strIdx+1;
    if i == 2
        strPlt{strIdx} = plotString('$N_e$',parSim.OP.engSpd*60,1,'rpm');
    elseif i == 3
        strPlt{strIdx} = plotString('$p_{Int}$',parSim.OP.pInt*1e-5,100,'bar');
    elseif i == 4
        strPlt{strIdx} = plotString('$p_{rail}$',parSim.OP.pRail*1e-5,1,'bar');
    elseif i == 5
        strPlt{strIdx} = plotString('$x_{bg}$',parSim.OP.xiBg*100,1,'\%');
    end
    
end

strIdx = strIdx+1;
strPlt{strIdx} = '';
strIdx = strIdx+1;
strPlt{strIdx} = '\textbf{Model Parameters, Combustion}';
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$w_{\rm min}$',parOpt.Reference.IMEP*1e-5,10,'bar');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$p_{\rm max}$',parOpt.Reference.pMax*1e-5,1,'bar');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$p''_{\rm max}$',parOpt.Reference.dpMax*1e-5,10,'${\rm bar}/^\circ \rm CA$');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$T_{\rm evo,min}$',parOpt.Reference.Tmin-273,1,'$^\circ \rm C$');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$\theta_{\rm coc,max}$',parOpt.Reference.CoCmax,1,'$^\circ \rm CA$');
if parOpt.enNOx
%     strIdx = strIdx+1;
%     strPlt{strIdx} = valS('$T_{\rm bz}$',parOpt.Reference.TbzMax,1,'K');
    strIdx = strIdx+1;
    strPlt{strIdx} = plotString('$X_{\rm NO, max}$',parOpt.Reference.cNOx,1,'ppm');
end

strIdx = strIdx+1;
strPlt{strIdx} = '';
strIdx = strIdx+1;
strPlt{strIdx} = '\textbf{Simulation Results}';
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$w_{hp}$',imep,10,'bar');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('${\rm max}(p)$',pMax,1,'bar');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('${\rm max}(p'')$',dpMax,10,'${\rm bar}/^\circ \rm CA$');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$T_{\rm evo}$',simOpt.y(1,end)-273,1,'$^\circ \rm C$');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$\theta_{\rm coc}$',CA50,10,'$^\circ \rm CA$');
if parOpt.enNOx
    strIdx = strIdx+1;
    strPlt{strIdx} = plotString('${\rm max}(T_{\rm bz})$',Tbz,1,'K');
    strIdx = strIdx+1;
    strPlt{strIdx} = plotString('$X_{\rm NO,out}$',NOx,1,'ppm');
end
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$\Phi$',Phi,100,'');


strIdx = strIdx+1;
strPlt{strIdx} = '';
strIdx = strIdx+1;
strPlt{strIdx} = '\textbf{Optimization Results}';
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$m_{\rm f,tot}$',simOpt.y(5,end)*1e6,1000,'mg$/$str');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$u_{\rm soe}$',uSOE,10,'$^\circ \rm CA$');
strIdx = strIdx+1;
strPlt{strIdx} = plotString('$u_{\rm doe}$',uDOE*1e6,1,'$\mu$s');
strIdx = strIdx+1;

hAnnotation = annotation('textbox',[.5 .05 .5 .93],'string',strPlt);
set(hAnnotation,'Interpreter','latex')


end


%% helper functions
function[roundVal] = rnd2(value,dec)
roundVal = round(value(:)'*dec)/dec;

end

function[strOut] = rndStr(value,fac)
strOut = num2str(rnd2(value,fac));

end

function[strOut] = plotString(name,value,fac,varargin)
if nargin == 4
    
    strOut = [name ' = ' rndStr(value,fac) ' '  varargin{1}];
    
else
    
    strOut = [name ' = ' rndStr(value,fac)];
    
end
end
