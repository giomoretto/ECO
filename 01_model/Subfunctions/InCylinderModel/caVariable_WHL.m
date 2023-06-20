function [dropFactor] = caVariable_WHL(ca,parModel,varargin)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

if nargin == 3
    plotFigure = varargin{1};
else
    plotFigure = false;
    
end

deltaCA = parModel.wall.dropPercentage/100/2;
ca01 = parModel.wall.caDropLocation;
ca02 = ca01+deltaCA;

% smoothness factor
epsSmooth = parModel.wall.smoothness;


dropFactor = 1 -(sqrt((ca-ca01).^2 + epsSmooth)) - (ca-ca01) + ...
    (sqrt((ca-ca02).^2+epsSmooth))+ (ca-ca02);


if plotFigure
    figure(1); hold on;
    plot(ca,dropFactor)
    xlim([ca(1) ca(end)])
end

end

