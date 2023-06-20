function [mAct,mtot,DOI] = AlgebraicInjectorModel(t,pRail,SOE,DOE)

  % Inputs:
    % t:        time vector [mus]
    % pRail:    Rail Pressure [Pa]
    % SOE:      Start of Injector Energyzing [mus]
    % DOE:      Duration of Injector Energyzing [mus]
  % Outputs
    % m_kg:     fuel mass trajectory for t vector [kg]

    
  % Parameter Definition for injection rate model
    rho_diesel   = 830;
    d_noz        = 120*10^-6;
    n_noz        = 8;
    A_noz_tot    = (d_noz/2)^2*pi*n_noz;
    inj_del      = 55; % 210;    % [mus]     SOE to SOI delay (hydraulic delay)
    DOEmin       = 66.73;     % [mus]
  % Limit DOE and convert to s
    epsLim = 1;
    DOElim = 1/2*sqrt((DOE-DOEmin).^2+epsLim) + 1/2*(DOE-DOEmin) + DOEmin;
    
  % Injection rate model
    % bernoulli speed
      P_Cyl        = 60e5; 
      delta_p      = pRail-P_Cyl;
      u_bernoulli  = sqrt(2*delta_p/rho_diesel);
    % injector opening & closing rate gradient [g/s / mus]
      % rail pressure dependent correction factor
        corr_fct_s  = 3.006-9090*pRail^-0.4791;
      % gradient at injector opening (g/s *1/mus)
        so =  100* corr_fct_s * 1e-3; % 80 * corr_fct_s * 1e-3;
      % gradient at injector closing (g/s *1/mus)
        sc =  200* corr_fct_s * 1e-3; % 130* corr_fct_s * 1e-3;
    % Real (physical) duration of injection
    % (start of opening to end of closing) [mus]
      DOI = (DOElim-DOEmin)*2.093; % (DOElim-DOEmin)*2.093;
    % discharge coefficient, dependent on DOI & rail pressure [-]
      corr_fct_cd = 9.8e-9*pRail+0.425;
      DOIlim      = 1/2*sqrt((DOI-70).^2+epsLim) + 1/2*(DOI-70) + 70; % mus
      coeff_d     = (0.9301 - 0.02744*(DOIlim*1e-6).^-0.3676)*corr_fct_cd;
      
    % limit discharge coefficient around doe = 100 mus
      DOIlower    = 70;
      eps         = 1e-2;
      coeff_dlim  = (0.9301 - 0.02744*(DOIlower*1e-6).^-0.3676)*corr_fct_cd;
      coeff_d     = 0.5*(coeff_d + coeff_dlim + sqrt(eps + (coeff_d-coeff_dlim).^2));
      
    % injection rate at fully open injector [kg/s] / [kg]
      % Injection flow rate according to bernoulli [g/s]
        dmFuel      = coeff_d.*u_bernoulli*A_noz_tot*rho_diesel*1e3;
      % limitation of flow rate for not completely open injectors
      % (trapezoid assumtion) [g/s]
        dmFuelMax   = DOI*(so*sc)/(so+sc);
        eps         = 1e-2;
        dmFuel      = dmFuelMax - ...
                     1/2*sqrt((dmFuel-dmFuelMax).^2 + eps) + ...
                     1/2*(dmFuel-dmFuelMax);
    % Complete fuel injection value [kg]
      cost_type   = 'exact'; % exact / quadratic / linear
      switch cost_type
          case 'exact'
              mtot  = sum(1e-9*(DOI.*dmFuel -1/2*dmFuel.^2*(1/so+1/sc)));
          case 'quadratic'
              mtot  = 1e-6*sum(0.02767*(DOE-100) + 9.004e-5*(DOE-100).^2);
          case 'linear'
              mtot  = 1e-6*sum(0.04286*(DOE-100));
          otherwise
              disp('use exact cost function')
              mtot  = sum(1e-9*(DOI.*dmFuel -1/2*dmFuel.^2*(1/so+1/sc)));
      end
    % Duration of Opening and Closing
      dt_start = dmFuel/so;
      dt_stop  = dmFuel/sc;
    % Derive SOI and EOI [mus]
      t_SOI    = SOE+inj_del;
      t_EOI    = t_SOI+DOI-dt_stop;
    % Derive smothness of function (Optimization only works with const. eps)
      eps11    = 500; 
      eps12    = 500; %   + (DOI*1e6).^2*1/40;
      eps21    = 500; % 1000;
      eps22    = 500; % 1000;
    % Fuel mass trajectory based on four "integrated" limit functions [kg]
      m    = 1e-9*(so/8)*(...
             (sqrt((t-t_SOI).^2+eps11)+(t-t_SOI)).^2 -...
             (sqrt((t-t_SOI-dt_start).^2+eps12)+(t-t_SOI-dt_start)).^2)...
             -1e-9*(sc/8)*(...
             (sqrt((t-t_EOI).^2+eps21)+(t-t_EOI)).^2 -...
             (sqrt((t-t_EOI-dt_stop).^2+eps22)+(t-t_EOI-dt_stop)).^2);
      % sum up all masses of individual injections
      mAct = sum(m);
end