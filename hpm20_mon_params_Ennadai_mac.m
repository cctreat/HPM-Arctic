% Parameter values for hpm20_mon;  S. Frolking, J. Talbot, 2008-2015

% v.20_mon: cleaning things up ? May 2018

% v.20: making version for version control - June 2015

% v.13: re-organized to make changing # of PFTs easier (Sept. 2014); set code to take arbitrary number
% of PFTs.  Added old/new Carbon (double PFTs).  

% v. 12: modified to read in atmosphere del-14C and more flexible simulation times; summer 2014

% v.8: adjusted lawn and hummock sphagnum NPP sensitivity to WTD

% new in version 6: WTD_range and PD_range have 2 sets of values (see 'vegetation_NPP_4.m')

% ********************************************************
% FILE BASE NAMEs AND SOME GENERAL SIMULATION FLAGS

% Build outfile name from model version, site, climate, years, other; add director pathways
%   have site and climate names used in climate file name.

site_name = 'Ennadai';
sim_name = '_1_NPP05';
monthly_T_P_name =  '_monthly_T_P_5810BP_2100CE'; 
working_directory = pwd;
dataWrite_workDirect = '~/Dropbox/Research/UNH Arctic HPM/Permafrost Gradient/Analysis/';

out_name = strcat(dataWrite_workDirect, 'hpm20_mon_output_files/', site_name, sim_name);
in_name = strcat(dataWrite_workDirect, 'hpm20_mon_input_files/', site_name, sim_name);
clim_in_name = strcat(dataWrite_workDirect, 'climate_drivers/',site_name, monthly_T_P_name,'.csv');
c14_in_name = strcat('~/Dropbox/HPM30_monthly_time_step/hpm20_mon_input_files/','annual_atm_del_14C_20000BP_to_2500AD_all_RCP','.csv');

sim_start = 5810; % years BP (before 'present'), where 0 BP = 1950 CE
sim_end = -150;   % years BP  (-150 BP = 2100 CE)
sim_len = sim_start - sim_end + 1;  % simulation length (years)

gipl_flag = 1; % if 0 (or 1) skip (or run) GIPL soil physics model: no (or yes) temperature effect on decomp
%   gipl_flag should always be 1??
RCP_flag = 1; % 1 = RCP8.5, 2 = RCP6.0, 3 = RDCP4.5, 4 = RCP2.6 (used for 21st century 14C values from Heather Graven)
pf_flag = 1; % if 1 site has or may sometimes have permafrost; otherwise 0 

max_pot_peat_ht = 6; % max. height for binning 'fancy' graphs

% **************
%  SITE CLIMATE

% read in monthly climate data from climate processing file
% hpm_climate_params20;
ann_temp = -9;  % site mean annual temp (C) for parameter values
ann_ppt = 0.29;  % site total annual precipitation (m/y) for parameter values
ann_ET_0 = 0.29;  % site base evapotranspiration (m/y) used to compute base run-off

latitude = 61; % degrees North > 0

% compute monthly mean daylength (fraction of 24 hours) using algorithm from WBM code
solstice	= 23.44/180 * pi;
month_midday = [15 46 74 105 135 166 196 227 258 288 319 349];  
dayLength = zeros(1,12);

for (imonth = 1:12)
    dec	= solstice * cos(2 * pi/365 * (month_midday(imonth)));
    arg	= -tan(dec) * tan(latitude/180*pi);
    arg	= -(arg<-1) + (arg>1) + (abs(arg)<=1)*arg;
    dayLength(imonth) = (1-acos(arg)/pi);
end

%  some initialization of threshold values

ald_0 = 1.0;  % first year active layer depth, if needed (m)
wtd_0 = 0.04; % initialization period water table depth (m)
start_depth = 0.25; % depth of initial peat accumulation (m) at which water balance calculations begin

% *********************
%  LATERAL HEAT FLUX PARAMETERS

HeatFlux_DeltaT = 4;  % degree C shift in thawed peat monthly temp (note 1°C into 1 m saturated peat is about 1.5 W/m2)
HeatFlux_StartYear = 250;  % first year of simulation with lateral heat flux 
HeatFlux_EndYear   = 1500;  % last year of simulation with lateral heat flux 

% ************************
%  DECOMP PARAMETERS

% initial decomposition (mass-loss) rates and anoxia factor 
%   (make anoxia factor more variable, as in new paper by Blodau?)

wfps_opt = 0.45;  % must be <= 0.5; optimum WFPS for decomposition (see spreadsheet 'simple decomp models.xls')
wfps_max_rate = 1.0;   % decomp rate multiplier at WFPS = WFPS_opt.
wfps_sat = 1.0;      % WFPS at saturation
wfps_sat_rate = 0.3;   % decomp rate multiplier at WFPS = 1.0 (i.e., at annual mean WTD).  Change to monthly (lower?)
wfps_min_rate = 0.001;   % decomp rate multiplier minimum, deep in catotelm.   Okay for monthly?
wfps_curve = (wfps_sat - wfps_opt)^2 / (4 * (wfps_max_rate - wfps_sat_rate)); % parabola w/ value = wfps_sat_rate at WFPS = 1.0


% **************
%  VEGETATION 

%   NOTE: plants don't ?grow? or accumulate biomass (change for trees [JT]?), so litterfall = NPP

% ARCTIC VERSION  
%    arctic version will use active layer depth rather than peat height for NPP

num_veg = 5;
 
PFT_names = [ 'moss' 'sedge_ag' 'sedge_bg' 'shrub_ag' 'shrub_bg']; 
 
mosses =    [ 1 0 0 0 0 ];
vasculars = [ 0 1 1 1 1 ];
sedges =    [ 0 1 1 0 0 ];
woody =     [ 0 0 0 1 1 ];

PFT_param = zeros(num_veg,12);
 
% *** PFT Parameters                   ** PD not used **
%                 WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
PFT_param(1,:) = [ 0.1   0.09    0.35   1.0   2.   19.   1.0    19.    29.     0.5      1.0     0.04  ]; % moss
PFT_param(2,:) = [ 0.025 0.15    0.20   1.0   2.   19.   1.5    1.0    29.     1.0      1.0     0.25  ]; % sedge aboveground
PFT_param(3,:) = [ 0.025 0.15    0.20   1.0   2.   19.   1.5    1.0    29.     1.0      0.0     0.225 ]; % sedge belowground
PFT_param(4,:) = [ 0.3   0.15    3.5    1.0   2.   19.   2.0    1.5    29.     1.3      1.0     0.15  ]; % shrub aboveground
PFT_param(5,:) = [ 0.3   0.15    3.5    1.0   2.   19.   2.0    1.5    29.     0.7      0.0     0.10  ]; % shrub belowground

%  'tf_xxx' parameters below are used if vascular PFTs are partitioned into two
%    component PFTs (aboveground and belowground), with some different parameters (e.g., k_exp)
%      generalize from sedge root to all roots?  Eliminate altogether?

% tf_sedge_root =     [ 0 0 0 1 0 ];
% tf_non_sedge_root = [ 0 0 0 0 1 ];
% tf_moss =           [ 1 0 0 0 0 ];
% tf_ag_pft =         [ 1 1 1 1 1 ] - tf_sedge_root - tf_non_sedge_root;

% NON-ARCTIC VERSION
%    non-arctic version will use peat height rather than active layer depth for NPP

% num_veg = 13;
% 
% PFT_names = ['min_grass' 'min_herb' 'min_sedge' 'decid_shrub' ?brown_moss' 'hollow_sphagnum' 'lawn_sphagnum' ...       ?                  'hummock_sphagnum' 'feather_moss' 'omb_herb' 'omb_sedge' 'evrgn_shrub' 'tree']; 
% 
% mosses =    [ 0 0 0 0 1 1 1 1 1 0 0 0 0 ];
% vasculars = [ 1 1 1 1 0 0 0 0 0 1 1 1 1 ];
% sedges =    [ 0 0 1 0 0 0 0 0 0 0 1 0 0 ];
% woody =     [ 0 0 0 1 0 0 0 0 0 0 0 1 1 ];

% % *** PFT Parameters                                    ***ALD NOT USED****
% %                  WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
% PFT_param(1,:) =  [ 0.4    0.4    0.4    0.01  1.0   1.0   0.1    1.     1.     3*0.75   0.5     0.2  ]; % minerotrophic grass
% PFT_param(2,:) =  [ 0.1    0.3    0.3    0.3   1.0   1.0   0.1    1.     1.     3*1.0    0.2     0.4  ]; % minerotrophic herb
% PFT_param(3,:) =  [ 0.1    0.4    0.4    0.1   2.0   2.0   0.1    1.     1.     3*1.0    0.2     0.3  ]; % minerotrophic sedge
% PFT_param(4,:) =  [ 0.2    0.2    1.0    1.0   2.0   2.0   0.1    1.     1.     3*0.5    0.5     0.25 ]; % deciduous shrub
% PFT_param(5,:) =  [ 0.01   0.2    0.05   0.1   1.5   1.5   0.1    1.     1.     1*0.5    1.0     0.1  ]; % brown moss
% PFT_param(6,:) =  [ 0.01   0.2    0.05   2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.1  ]; % hollow sphagnum
% PFT_param(7,:) =  [ 0.1    0.3    0.4    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.07 ]; % lawn sphagnum
% PFT_param(8,:) =  [ 0.2    0.1    0.5    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.05 ]; % hummock sphagnum
% PFT_param(9,:) =  [ 0.4    0.4    0.6    4.0   6.0   19.   0.1    1.     1.     1*0.25   1.0     0.1  ]; % feathermoss
% PFT_param(10,:) = [ 0.2    0.2    0.2    4.0   2.0   19.   0.1    1.     1.     1*0.25   0.5     0.3  ]; % ombrotrophic herb
% PFT_param(11,:) = [ 0.2    0.3    0.3    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.2     0.2  ]; % ombrotrophic sedge
% PFT_param(12,:) = [ 0.3    0.3    1.0    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.5     0.2  ]; % evergreen shrub
% PFT_param(13,:) = [ 0.8    0.3    10.0   2.0   20.   10.   0.1    1.     1.     1*2.8    0.7     0.3  ]; % tree

% REDUCING NON-ARCTIC PFT NUMBER TO TRY OLD/NEW C RUN AT MB

% num_veg = 10;
% 
% PFT_1_name = ['min_sedge' 'decid_shrub' ?hollow_sphagnum? 'lawn_sphagnum? 'hummock_sphagnum' 'feather_moss' ...
%                  'omb_herb' 'omb_sedge' 'evrgn_shrub' 'tree']; 
% 
% mosses =    [ 0 0 1 1 1 1 0 0 0 0 ];
% vasculars = [ 1 1 0 0 0 0 1 1 1 1 ];
% sedges =    [ 1 0 0 0 0 0 0 1 0 0 ];
% woody =     [ 0 1 0 0 0 0 0 0 1 1 ];

% % *** PFT Parameters                                    ***ALD NOT USED****
% %                  WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
% PFT_param(1,:) =  [ 0.1    0.4    0.4    0.1   2.0   2.0   0.1    1.     1.     3*1.0    0.2     0.3  ]; % minerotrophic sedge
% PFT_param(2,:) =  [ 0.2    0.2    1.0    1.0   2.0   2.0   0.1    1.     1.     3*0.5    0.5     0.25 ]; % deciduous shrub
% PFT_param(3,:) =  [ 0.01   0.2    0.05   2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.1  ]; % hollow sphagnum
% PFT_param(4,:) =  [ 0.1    0.3    0.4    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.07 ]; % lawn sphagnum
% PFT_param(5,:) =  [ 0.2    0.1    0.5    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.05 ]; % hummock sphagnum
% PFT_param(6,:) =  [ 0.4    0.4    0.6    4.0   6.0   19.   0.1    1.     1.     1*0.25   1.0     0.1  ]; % feathermoss
% PFT_param(7,:) = [ 0.2    0.2    0.2    4.0   2.0   19.   0.1    1.     1.     1*0.25   0.5     0.3  ];  % ombrotrophic herb
% PFT_param(8,:) = [ 0.2    0.3    0.3    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.2     0.2  ];  % ombrotrophic sedge
% PFT_param(9,:) = [ 0.3    0.3    1.0    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.5     0.2  ];  % evergreen shrub
% PFT_param(10,:) = [ 0.8    0.3    10.0   2.0   20.   10.   0.1    1.     1.     1*2.8    0.7     0.3  ]; % tree

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
% DECOMPOSITION- LITTER QUALITY  

% **************
% ARCTIC VERSION  

k_exp_temp = 4.;   % litter incubation temperature (°C) from Hobbie paper

k_0 = k_exp .* (1 + 3 * k_exp);  %see spreadsheet 'simple decomp models.xls'; adjusts k_0 for m/m0 model of decay
k_month_0 = k_0 / 0.411 / 12;  % convert 'per-0.4-yr' (Hobbie study) to 'per-year' to 'per-month'

% ???   does this make sense (next line)?????
% k_month_0 = k_0 / 12 * (2 + 3*exp(-ann_temp/9)).^(-(ann_temp - 10)/10);  % convert year at 6° (mer bleue) to year at 10° (ref)
k_month_0 = k_month_0 * (2 + 3*exp(-k_exp_temp/9)).^(-(k_exp_temp - 10)/10);  % convert month at 4° (Hobbie incubation) to year at 10° (ref)


% **************
% NON-ARCTIC VERSION  

% tt = 5; % years of decomp to match algorithms (5 is reasonable, i.e., long for litter bags)
% k_0 = k_exp * (1 + ((exp(k_exp*tt) - 1) ./ (k_exp * tt) -1);
% k_0_month = k_0 / 12;  % will this work okay?
% k_month_0 = k_0 / 12 * (2 + 3*exp(-ann_temp/9)).^(-(ann_temp - 10)/10);  % convert year at 6° (mer bleue) to year at 10° (ref)

%--- BUILD TOTAL NPP SURFACE
%     No Permafrost (pf_flag = 0) uses WTD and PD
%      Permafrost (pf_flag = 1) uses WTD and ALD

% ONLY MOSSES  (for sensitivity testing)
% NPP_rel   = NPP_rel .* mosses;

% ONLY VASCULAR  (for sensitivity testing)
% NPP_rel   = NPP_rel .* vascular;

if (pf_flag < 0.5)
    total_npp = hpm20_mon_vegNPP(WTD_opt,WTD_range,PD_opt,PD_range,NPP_rel,k_0,num_veg); % using peat depth
else
    total_npp = hpm20_mon_vegNPP(WTD_opt,WTD_range,ALD_opt,ALD_range,NPP_rel,k_0,num_veg); % using active layer depth
end

% Specify site absolute maximum NPP (kg/m2/y dry matter) during peatland lifetime

max_npp = 0.5;   % approximate absolute maximum total NPP for all vegetation at mean annual T = 10°C, kg/m2/y
                         %  for TOOLIK (ann_temp = -10°C) Q10 multiplier is 1.5^(-2) = 0.44
% original value was 1
q10_npp = 1.5;   % see Julie Talbot email of 4 June 2014
max_npp = max_npp * q10_npp^((ann_temp - 10)/10);

NPP_rel = NPP_rel * (max_npp / total_npp)   % scale relative NPP of all PFTs so that max sum NPP ~ 'max_npp'

% # years averaging WTD for vascular plant NPP (1 year for non-vascular)
%   ?? add another lag value for trees different from other vascular?

lag_years = 10;  


% ********************************************************
% RUN WITH DOUBLE PFTS FOR OLD-NEW CARBON ANALYSIS

tf_old_new = 0; % 1: double PFTs for old/new; otherwise = 0 & do not do this
tf_old_new_timing = 150;  % years before end of simulation to switch 

if (tf_old_new > 0.5)
    year_old2new = sim_len - tf_old_new_timing;
    num_veg = 2. * num_veg;
    PD_opt = [PD_opt PD_opt];
    PD_range = [PD_range PD_range];
    ALD_opt = [ALD_opt ALD_opt];
    ALD_range = [ALD_range ALD_range];
    WTD_opt = [WTD_opt WTD_opt];
    WTD_range = [WTD_range WTD_range];
    ag_frac_npp = [ag_frac_npp ag_frac_npp];
    bg_frac_npp = [bg_frac_npp bg_frac_npp];
    k_0 = [k_0 k_0];
    k_month_0 = [k_month_0 k_month_0];
    NPP_rel = [NPP_rel NPP_rel];
    mosses =    [ mosses mosses ];
    vasculars = [ vasculars vasculars ];
    sedges =    [ sedges sedges ];
    woody =     [ woody woody ];
%     tf_sedge_root = [tf_sedge_root tf_sedge_root];
%     tf_non_sedge_root = [tf_non_sedge_root tf_non_sedge_root];
%     tf_moss = [tf_moss tf_moss];
%     tf_ag_pft = [tf_ag_pft tf_ag_pft];
else
    year_old2new = sim_len + tf_old_new_timing;  % i.e., never happens
end


% ****************************
% ROOT PROFILES AND LITTER INPUT

% FOLLOWING BAUER (2004)
%   sedge root litter input: exponential decay profile to 2 m; 80% above 0.2 to 0.3 m.
%          litter_input [kg/m2/m] = beta * exp(-alpha * depth) * layer_thickness
%   other vascular root litter: ~constant profile to max(WTD, 0.2 m)
%          litter_input [kg/m2/m] = [bg_npp / max(WTD, 0.2)] * layer_thickness
%   moss root litter input equals zero

rootin_min = 0.2;  % minimum depth [m] for non-sedge root profile (most roots above max of WTD & rootin_min)
rootin_max = 1.0; % max depth (m) for non-sedge roots
rootin_d80 = 0.3;   % depth in meters to 80% of sedge roots       
rootin_alpha = -log(rootin_d80) / (1. - 0.8); % controls sedge root profile expon. decay; denom. = 1 - frac. roots above 'd80'
rootin_sedge_max = 2.0;   % maximum depth [m] for sedge root profile
rootin_c5 = 0.04;    % no longer used, this was a smoothing term for the root distribution


% **************
%   BULK DENSITY

min_bulk_dens = 50.;   % kg/m3
del_bulk_dens = 80.;   % bulk density increase down profile, kg/m3
dens_c1 = 0.333;  % m_star value at which bulk density rises halfway from min to max
dens_c2 = 0.20;  % parameter controlling steepness of bulk density transition (smaller is steeper)
OM_dens = 1300; % density of organic matter [kg/m3]


% ****************************
%  WATER BALANCE

del_water_threshold = 0.0002; % MINIMUM change in peat water (m) to recompute monthly water profile

% **************

%  BOG or FEN or PERMAFROST ??

bog_fen_id = 3;   % 1= fen-to-bog. 2 = persistent fen, 3 = permafrost

if (bog_fen_id < 1.5)   % FEN-TO-BOG VALUES
    
    Roff_c1 = max(0,ann_ppt - ann_ET_0 + 0.1); % + 0.05 % max runoff + max ET = mean annual precip + 0.05 m/y 
    Roff_c2 = 0.2;  % linear increase in runoff (m/y) per meter of total peat height
    Roff_c2a = 1.2 * start_depth;  % peat height needed to get base run-off (factor = 1 + c2*(H-c2a))

    anoxia_scale_length = 0.3;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate

    runon_c1 = 1.0;  % total peat height (m) at which run-on declines by ~50%
    runon_c2 = 0.5;  % controls rate of decline of run-on as function of peat height (see 'HPM vegetation productivity.xls')
    runon_c3 = 0.0;   % magnitude of maximum run-on (m/month)

elseif (bog_fen_id < 2.5)   %  PERENNIAL FEN VALUES
    
    Roff_c1 = max(0,ann_ppt - ann_ET_0 + 0.1); % max runoff + max ET = mean annual precip
    Roff_c2 = 0.02;  % linear increase in runoff (m/y) per meter of total peat height
    Roff_c2a = 1.;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))

    anoxia_scale_length = 3.0;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate

    runon_c1 = 5.0;  % total peat depth (m) at which run-on declines by ~50%
    runon_c2 = 0.5;  % controls rate of decline of run-on as function of peat height (see 'HPM vegetation productivity.xls')
    runon_c3 = 0.05; % magnitude of maximum run-on (m/month)

else   %  PERMAFROST SITE VALUES
    
    Roff_c1 = max(0,ann_ppt - ann_ET_0 + 0.25); % max runoff + max ET = mean annual precip + xx m/yr
    Roff_c2 = 0.2;  % linear increase in runoff (m/y) per meter of total peat height
    Roff_c2a = 1.;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))
    Roff_c2a = 1.2 * start_depth;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))

    anoxia_scale_length = 3;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate

    runon_c1 = 0.5;  % total peat depth (m) at which run-on declines by ~50%
    runon_c2 = 0.5;  % controls rate of decline of run-on as function of peat height (see 'HPM vegetation productivity.xls')
    runon_c3 = 0.05;  % magnitude of maximum run-on (m/month)

end
% ???  ***SHOULD Roff_c1 EVER BE ZERO???***

% for value of Roff_c1 see top of this script.
Roff_c3 = 0.25; % minimum profile relative transmissivity (see hpm_WatBal7.m')
Roff_c4 = -0.1; % threshold water table depth for extra spillover (= Roff_c4 - WTD)
max_inundation = 0.02;  % 0.1; % max. water pooling depth (m) ABOVE peat surface

% modifications June 2008 based on Lafleur et al. 2005 ET from Mer Bleue paper
ET_wtd_1 = 0.15;   % WTD threshold for full ET (m)
ET_wtd_2 = 0.5;   % WTD threshold for low ET (m)
ET_min_frac = 0.2;   % min ET factor (fraction of PET)
ET_param = 1/(ET_wtd_2 - ET_wtd_1);   % linear drop in ET as WTD drops from ET_wtd_1 to ET_wtd_2
ET_snow_depth = 0.05; % monthly mean snowdepth (m) for which AET goes to zero (snow sublimation?

% Parameters (see excel spreadsheet: 'anoxia & bulk dens & WFPS % profile.xls')
% calculates peat cohort fractional water content (0 - 1) above the water table
%   depends on distance above water table and peat bulk density;

wfps_c1 = 0.03;
wfps_c2 = 0.5;
wfps_c3 = 20;


% ****************************
%  RADIOCARBON

tau_c14 = 5730. / log(2); % radiocarbon decay rate per year ('log' in matlab is natural logarithm)


% write a few things to the terminal  (keep this?)
disp(sprintf('numveg: %d  total NPP (kg/m2/y): %d', num_veg, total_npp));
disp(sprintf('ann_ppt (m/y): %d  ann_ET0 (m/y): %d, ann_runoff0 (m/y): %d', ann_ppt, ann_ET_0, Roff_c1));

% ****************************
%  SAVE PARAMETERS

save('hpm20_mon_param_vals','out_name', 'in_name', 'clim_in_name', 'c14_in_name', 'site_name', 'sim_name', 'monthly_T_P_name', ...
    'sim_len','sim_start','sim_end','tau_c14','max_pot_peat_ht',...
    'gipl_flag','pf_flag','RCP_flag', 'latitude', 'dayLength', 'month_midday', ...
    'start_depth', 'ald_0', 'wtd_0', 'lag_years', ...
    'num_veg','mosses','vasculars','sedges','woody', 'tf_old_new','year_old2new', ...
    'WTD_opt','WTD_range','PD_opt','PD_range','ALD_opt','ALD_range',...
    'NPP_rel','max_npp', 'ag_frac_npp','bg_frac_npp','q10_npp', ...
    'rootin_d80','rootin_alpha','rootin_min','rootin_sedge_max','rootin_c5','rootin_max',...
    'k_0','k_month_0','wfps_opt','wfps_curve','wfps_sat_rate','wfps_min_rate','anoxia_scale_length',...
    'min_bulk_dens','del_bulk_dens','dens_c1','dens_c2','OM_dens',...
    'ann_temp','ann_ppt', 'ann_ET_0','del_water_threshold',...
    'ET_wtd_1','ET_wtd_2','ET_min_frac','ET_param','ET_snow_depth',...
    'Roff_c1','Roff_c2','Roff_c2a','Roff_c3','Roff_c4','max_inundation','runon_c1','runon_c2','runon_c3',...
    'wfps_c1','wfps_c2','wfps_c3','HeatFlux_DeltaT', 'HeatFlux_StartYear','HeatFlux_EndYear');

