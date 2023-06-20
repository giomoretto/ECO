function [kappaTot] = calcKappa(xi,theta,specR)

coef = struct();
% Table data with polynomial coefficients for each species
% Source: http://combustion.berkeley.edu/gri-mech/version30/files30/thermo30.dat
% 1:7  --> a1 - a7 for theta  > 1000 (high temperature regime)
% 8:14 --> a1 - a7 for theta <= 1000 (low temperature regime)
% N2
coef.N2 =...
    [ 0.02926640E+02;...
      0.14879768E-02;...
     -0.05684760E-05;...
      0.10097038E-09;...
     -0.06753351E-13;...
     -0.09227977E+04;...
      0.05980528E+02;...
      0.03298677E+02;...
      0.14082404E-02;...
     -0.03963222E-04;...
      0.05641515E-07;...
     -0.02444854E-10;...
     -0.10208999E+04;...
      0.03950372E+02];

% O2
coef.O2 =...
    [ 3.28253784E+00;...
      1.48308754E-03;...
     -7.57966669E-07;...
      2.09470555E-10;...
     -2.16717794E-14;...
     -1.08845772E+03;...
      5.45323129E+00;...
      3.78245636E+00;...
     -2.99673416E-03;...
      9.84730201E-06;...
     -9.68129509E-09;...
      3.24372837E-12;...
     -1.06394356E+03;...
      3.65767573E+00];
  
% CO2
coef.CO2 =...
    [ 3.85746029E+00;...
      4.41437026E-03;...
     -2.21481404E-06;...
      5.23490188E-10;...
     -4.72084164E-14;...
     -4.87591660E+04;...
      2.27163806E+00;...
      2.35677352E+00;...
      8.98459677E-03;...
     -7.12356269E-06;...
      2.45919022E-09;...
     -1.43699548E-13;...
     -4.83719697E+04;...
      9.90105222E+00];
  
% H2O
coef.H2O =...
    [ 3.03399249E+00;...
      2.17691804E-03;...
     -1.64072518E-07;...
     -9.70419870E-11;...
      1.68200992E-14;...
     -3.00042971E+04;...
      4.96677010E+00;...
      4.19864056E+00;...
     -2.03643410E-03;...
      6.52040211E-06;...
     -5.48797062E-09;...
      1.77197817E-12;...
     -3.02937267E+04;...
     -8.49032208E-01];

% C7H16, n-Hepton --> source: lect. notes 'Technical Combustion' ITV
coef.CxHy =...
    [ 0.22818893E+02;...
      0.32543454E-01;...
     -0.11120041E-04;...
      0.17131743E-08;...
     -0.96212101E-13;...
     -0.33678738E+05;...
     -0.94335007E+02;...
      0.30149546E+01;...
      0.54457203E-01;...
      0.21812681E-04;...
     -0.54234111E-07;...
      0.20808730E-10;...
     -0.26003379E+05;...
      0.17508575E+02];
  
% C8H18, iso-Octane
% coef.CxHy =...
%     [ 2.71373590E+01;...
%       3.79004890E-02;...
%      -1.29437358E-05;...
%       2.00760372E-09;...
%      -1.16400580E-13;...
%      -4.07958177E+04;...
%      -1.23277495E+02;...
%      -4.20868893E+00;...
%       1.11440581E-01;...
%      -7.91346582E-05;...
%       2.92406242E-08;...
%      -4.43743191E-12;...
%      -2.99446875E+04;...
%       4.49521701E+01];
     
species = fieldnames(xi);

cpTot    = 0;
specRTot = 0;
for k = 1:length(species)
    if theta > 1000
        aK = coef.(species{k})(1:5);
    else
        aK = coef.(species{k})(8:12);
    end
    
    cpTot = cpTot +...
        (aK(1) + aK(2)*theta + aK(3)*theta^2 + aK(4)*theta^3 + aK(5)*theta^4)*...
        specR.(species{k})*xi.(species{k});
    
    specRTot = specRTot + specR.(species{k})*xi.(species{k});
end

cvTot = cpTot - specRTot;
kappaTot = cpTot/cvTot;

end