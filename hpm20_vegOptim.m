function [WTD_opt, WTD_range, ALD_opt, ALD_range,...
    NPP_rel, max_npp, q10_npp, ag_frac_npp, bg_frac_npp, k_0,k_month_0] ...
    = hpm20_vegOptim(x, params)
% optimization function for parameters

% parameters to be optimized:
% PFTs: NPP relative
% height fen/bog transition and where run on turns to runoff

% **************
%  VEGETATION 

%   NOTE: plants don't ?grow? or accumulate biomass (change for trees [JT]?), so litterfall = NPP

% ARCTIC VERSION  
%    arctic version will use active layer depth rather than peat height for NPP

num_veg = params.num_veg;
 
PFT_param = zeros(num_veg,12);
 
% *** PFT Parameters                   ** PD not used **
%                 WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
PFT_param(1,:) = [ 0.1   0.09    0.20    1.0   2.   19.   1.0    19.    29.     1.0      1.0     0.04  ]; % moss
PFT_param(2,:) = [ 0.025 0.05    0.15    1.0   2.   19.   0.5    1.0    2.0     2.0      1.0     0.25  ]; % sedge aboveground
PFT_param(3,:) = [ 0.025 0.05    0.15    1.0   2.   19.   0.5    1.0    2.0     2.0      0.0     0.225 ]; % sedge belowground
PFT_param(4,:) = [ 0.25   0.09    1.5     1.0   2.   19.   0.8    1.0    9.     1.3*x(3)      1.0     0.15  ]; % shrub aboveground
PFT_param(5,:) = [ 0.25   0.09    1.5     1.0   2.   19.   0.8    1.0    9.     0.7*x(3)      0.0     0.10  ]; % shrub belowground



% **************
% VEGETATION ALL VERSIONS

for ii = 1:1:num_veg
    WTD_opt(ii)     = PFT_param(ii,1);      % optimal water table depth (m)
    WTD_range(1,ii) = PFT_param(ii,2);      % variance on shallow WTD side (m)
    WTD_range(2,ii) = PFT_param(ii,3);      % variance on deep WTD side (m)
    PD_opt(ii)      = PFT_param(ii,4);      % optimal peat height (m)
    PD_range(1,ii)  = PFT_param(ii,5);      % variance on shallow PD side (m)
    PD_range(2,ii)  = PFT_param(ii,6);      % variance on deep PD side (m)
    ALD_opt(ii)     = PFT_param(ii,7);      % optimal active layer depth (m)
    ALD_range(1,ii) = PFT_param(ii,8);      % variance on shallow ALD side (m)
    ALD_range(2,ii) = PFT_param(ii,9);      % variance on deep ALD side (m)
    NPP_rel(ii)     = PFT_param(ii,10);     % relative NPP
    ag_frac_npp(ii) = PFT_param(ii,11);     % aboveground fraction of NPP
    bg_frac_npp(ii) = 1. - ag_frac_npp(ii); % belowground fraction of NPP
    k_exp(ii)       = PFT_param(ii,12);     % litterbag k-value (1/year)
end


% **************
% ARCTIC VERSION  

k_exp_temp = params.k_exp_temp; %4.;   % litter incubation temperature (°C) from Hobbie paper

k_0 = k_exp .* (1 + 3 * k_exp);  %see spreadsheet 'simple decomp models.xls'; adjusts k_0 for m/m0 model of decay
k_month_0 = k_0 / 0.411 / 12;  % convert 'per-0.4-yr' (Hobbie study) to 'per-year' to 'per-month'
% k_0 = (k_exp.*1.0633);

% ???   does this make sense (next line)?????
% k_month_0 = k_0 / 12 * (2 + 3*exp(-ann_temp/9)).^(-(ann_temp - 10)/10);  % convert year at 6° (mer bleue) to year at 10° (ref)
k_month_0 = k_month_0 * (2 + 3*exp(-k_exp_temp/9)).^(-(k_exp_temp - 10)/10);  % convert month at 4° (Hobbie incubation) to year at 10° (ref)

%--- BUILD TOTAL NPP SURFACE
%     No Permafrost (pf_flag = 0) uses WTD and PD
%      Permafrost (pf_flag = 1) uses WTD and ALD

% ONLY MOSSES  (for sensitivity testing)
% NPP_rel   = NPP_rel .* mosses;

% ONLY VASCULAR  (for sensitivity testing)
% NPP_rel   = NPP_rel .* vascular;

total_npp = hpm20_mon_vegNPP(WTD_opt,WTD_range,ALD_opt,ALD_range,NPP_rel,k_0,num_veg); % using active layer depth

% Specify site absolute maximum NPP (kg/m2/y dry matter) during peatland lifetime

max_npp = 1.1;   % approximate absolute maximum total NPP for all vegetation at mean annual T = 10°C, kg/m2/y
                         %  for TOOLIK (ann_temp = -10°C) Q10 multiplier is 1.5^(-2) = 0.44
% original value was 1
q10_npp = params.q10_npp;   % see Julie Talbot email of 4 June 2014
max_npp = max_npp * q10_npp^((params.ann_temp - 10)/10);
NPP_rel = NPP_rel * (max_npp / total_npp);   % scale relative NPP of all PFTs so that max sum NPP ~ 'max_npp'

end




