function [mon_AET, mon_PET] = hpm20_mon_aet(WTD, SnowDepth, layer_frozen, monthly_air_temp, imonth, num_days_in_month, params)

% function calculates monthly ET (m); 
% note: snow melt and snow sublimation are computed separately, so set to zero if snowpack > threshold
% follows revised version of PAM from Roulet
% reference for PAM is:
% Hilbert D, NT Roulet, TR Moore. 2000. Modelling and analysis of peatlands as dynamical systems, J. Ecology. 88:230-242.

% ------------

% first compute monthly PET follow Hamon function

% compute monthly PET using Hamon Function (taken from WBM code provided by Danielle Grogan)
%   ASK DANIELLE AND DOMINIK ABOUT SNOW SUBLIMATION AND HAMON PET IN COLD WEATHER?
%   Reduce snow sublimation (now 10% per month?) by a daylength factor??
       
P_sat = (monthly_air_temp >= 0) * 0.61078*exp(17.26939 * monthly_air_temp / (monthly_air_temp + 237.3)) + ...
                   (monthly_air_temp < 0) * 0.61078*exp(21.874565 * monthly_air_temp / (monthly_air_temp + 265.5));
rho_sat = 2.167 * P_sat / (monthly_air_temp + 273.15);
PET_hamon = 0.3302 * params.dayLength(imonth) * rho_sat * num_days_in_month;   % m/month
mon_PET = PET_hamon;
     
% ------------

% AET PARAMETERS (values read in params script)
%      modify ALL for monthly?
%  ? ET_wtd_1: WTD threshold (m) for AET = PET
%  ? ET_wtd_2: WTD threshold (m) for AET = minimum
%  ? ET_min: minimum AET (m/month)
%  ? ET_param: used in annual model

if (WTD <= params.ET_wtd_1)
    ET = PET_hamon;
elseif (WTD >= params.ET_wtd_2)
    ET = params.ET_min_frac * PET_hamon;
else
    ET = PET_hamon * (params.ET_min_frac + (1 - params.ET_min_frac) * (WTD - params.ET_wtd_2) / (params.ET_wtd_1 - params.ET_wtd_2));
end

mon_AET = ET;

if (SnowDepth > params.ET_snow_depth)
    mon_AET = 0;
end

if (sum(layer_frozen(1:3)) >= 2)
    mon_AET = 0;
end

% END AET FUNCTION

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

