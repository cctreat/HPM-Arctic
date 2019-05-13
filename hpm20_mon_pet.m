function [mon_PET] = hpm20_mon_pet(monthly_air_temp, imonth, num_days_in_month, params)

% function calculates monthly HAMON PET (m); 
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

% END AET FUNCTION

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

