function rmseErr = hpm20_mon_mainOptim(x) %start runon/runoff = x(1), ...
% height_fbt = x(2); NPP_max = x(3); 
% anoxia_scalelengthFen = x(4); anoxia_scalelengthBog = x(5);


% HPM20_mon.m
% Monthly time-step version of HPM in matlab, based on HPM v.20
% Steve Frolking, May 2018

% functions called:  
%   hpm20_mon_params_*
%* VEGETATION AND PEAT
%   hpm20_mon_vegNPP  (called from hpm_params20 to initialize NPP)
%   hpm20_mon_npp
%   hpm20_mon_dens
%   hpm20_mon_decomp
%   hpm20_mon_rootin
%* WATER
%   hpm20_mon_wtd_wfps
%   hpm20_mon_wfps_only
%   hpm20_mon_runoff
%   hpm20_mon_runon
%   hpm20_mon_aet
%* SNOW AND ICE
%   hpm20_mon_snowpack
%   hpm20_mon_activelayer
%* SOIL TEMPERATURE
%   hpm20_gipl2_daily    UAF GI permafrost and soil temperature model from Marchenko
%   hpm20_gipl_params_pf          parameters for UAF GI model when simulation has permafrost (max depth 120 m)
%   hpm20_gipl_params_no_pf       parameters for UAF GI model when simulation doesn't have permafrost (max depth 10 m)
%   hpm20_mon_activelayer2        calculates active layer thickness (< 0 C for more than 2 years)
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% READ IN PARAMETERS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% hpm20_mon_params_Toolik;
% hpm20_mon_params_Lakkasuo;
% hpm20_mon_params_Seida;
%  hpm20_mon_params_Ennadai_win;
% hpm20_mon_params_Ennadai_mac
% hpm20_mon_params_JoeyL_winOptim;
% hpm20_mon_params_CCRP_win;
% hpm20_mon_params_Kukjuk_mac;
params=load('hpm20_mon_param_vals');
params.Roff_c2a = x(1);
params.runon_c1 = x(1);
params.anoxia_scale_length = x(4);

nveg = params.num_veg;
sim_len_yr = params.sim_len

if (params.tf_old_new > 0.5)  % using Old-New carbon switch
    old_new_ones = ones(1, nveg/2);
    old_new_zeros = zeros(1,nveg/2);
end

%--------------------------------------------
%  move GIPL outside of routine to speed up the model
%--------------------------------------------
% if (params.pf_flag > 0.5) % simulation with permafrost
%     hpm20_mon_gipl_params_pf;   % parameters for UAF GIPL 2.0 model with permafrost
% else
%     hpm20_mon_gipl_params_no_pf;   % parameters for UAF GIPL 2.0 model without permafrost
% end
% 
params_gipl = load('hpm20_mon_gipl_param_vals');   

ST1 = load('Monthly_soil_node_temp_matrix'); %have to be a little careful here, just uses soilT from last run.
soil_node_temp_month_mat = ST1.soil_node_temp_month_saveMAT;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% INITIALIZE ARRAYS   (*** use a consistent naming format, e.g., ?ann_name?, ?mon_name? ***)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% pre-allocate arrays to speed up simulations

% arrays by cohort and veg type

% small m arrays are masses as annual cohort (or layer) by veg types
m = zeros(sim_len_yr ,nveg);   % remaining mass in cohort (layer) i and veg type
m_0 = zeros(sim_len_yr ,nveg); % total input mass in cohort i and veg type
m_star = zeros(sim_len_yr ,nveg);   % = m / m_0
c14_m = zeros(sim_len_yr ,nveg);  % 14-C in each PFT of each cohort
k_mon = zeros(sim_len_yr ,nveg);  % mass loss rate (1/y)
rootin = zeros(sim_len_yr ,nveg);
ann_resp_array = zeros(sim_len_yr ,nveg);

% capital M vectors are masses as annual cohort accumulated across the veg types
M = zeros(sim_len_yr ,1);        % = sum m across veg. types in cohort i
M_0 = zeros(sim_len_yr ,1);      % = sum m_0 across veg. types in cohort i
M_star = zeros(sim_len_yr ,1);        % = M / M_0 in cohort i
M_overlying = zeros(sim_len_yr ,1);   % = sum M_total in profile above cohort/layer i
ann_del_M_tot = zeros(sim_len_yr ,1);   % annual change in total peat mass
c14_M = zeros(sim_len_yr ,1);    % 14-C in each cohort

% vectors down the annual peat cohort profile

depth = zeros(sim_len_yr ,1);  % cohort (layer) midpoint depth below peat surface in meters
thick = zeros(sim_len_yr ,1);  % cohort (layer) thickness in meters
zbottom = zeros(sim_len_yr ,1); % depth (m) from top of peat to bottom of cohort
porosity = zeros(sim_len_yr ,1);  % cohort (layer) porosity (m3/m3)
prev_thick = zeros(sim_len_yr ,1);  % cohort (layer) thickness in meters (from previous time step)
dens = zeros(sim_len_yr ,1);   % cohort (layer) bulk density in kg/m3
k_mean = zeros(sim_len_yr ,1);        % mass-weighted mean decomposition factor by cohort
anoxiafact = zeros(sim_len_yr ,1);    % anoxia profile, function of water table depth (anything else?)
wfps = zeros(sim_len_yr ,1);   % change to wfps

% these are temporary arrays, mostly used for shifting the peat profile down each year

m_temp = zeros(sim_len_yr ,nveg);
m_0_temp = zeros(sim_len_yr ,nveg);
m_star_temp = zeros(sim_len_yr ,nveg);
k_temp = zeros(sim_len_yr ,nveg);   % is this needed
rootin2 = zeros(sim_len_yr ,nveg);  % temporary array  % ## change to rootintemp
% agebiastemp = zeros(sim_len_yr ,nveg);
% agebiastemparr = zeros(sim_len_yr ,nveg);
depth_temp = depth;  % ## change to depthtemp
dens_temp = dens;   % ## change to denstemp
dens_prev = dens;    % used to make sure bulk density only increases over time (correct?) 
dens_prev_temp = dens_prev;   % ## change to dens_prevtemp

% arrays by time and veg type

ann_npp = zeros(sim_len_yr ,nveg);   % matrix of annual NPP by PFT (biomass units; in carbon units = biomass/2) 

% vectors by time 

time = zeros(sim_len_yr ,1); 
time_month = zeros(sim_len_yr * 12 ,1);
ann_NPP = zeros(sim_len_yr ,1);     % total annual NPP summed across PFTs (in biomass units; in carbon units = biomass/2) 
ann_RESP_C = zeros(sim_len_yr ,1);  % annual mass loss (in carbon units = biomass/2)   % ## change to ann_RESP_C
ann_RESP_C_old = zeros(sim_len_yr ,1);  % annual mass loss (in carbon units = biomass/2)   % ## change to ann_RESP_C
ann_RESP_C_new = zeros(sim_len_yr ,1);  % annual mass loss (in carbon units = biomass/2)   % ## change to ann_RESP_C
ann_ROOTIN = zeros(sim_len_yr ,1);      % ## change to ann_ROOTIN
ann_ROOTNPP = zeros(sim_len_yr ,1);   % ## change to ann_ROOTNPP
ann_AGMASSIN = zeros(sim_len_yr ,1);   % ## change to ann_AGMASSIN
ann_Z_total = zeros(sim_len_yr ,1);   % ## change to ann_Z_total
ann_del_Z_total = zeros(sim_len_yr ,1);   % ## change to ann_del_Z_total
ann_M_total = zeros(sim_len_yr ,1);   % ## change to ann_M_total
ann_RESP_del_c14 = zeros(sim_len_yr ,1); % del_14c of annual respiration   % ## change to ann_resp_del_C14

ann_temp_forcing = zeros(sim_len_yr ,1);       % deg C
ann_precip_forcing = zeros(sim_len_yr ,1);   % meters/year
ann_snowfall = zeros(sim_len_yr ,1);   % meters/year
ann_rainfall = zeros(sim_len_yr ,1);   % meters/year
ann_snowmelt = zeros(sim_len_yr ,1);   % meters/year
ann_snowsublimation = zeros(sim_len_yr ,1);   % meters/year
ann_max_snow_depth = zeros(sim_len_yr ,1);   % meters
ann_aet = zeros(sim_len_yr ,1);   % meters/year
ann_pet = zeros(sim_len_yr ,1);   % meters/year
ann_runon = zeros(sim_len_yr ,1);   % meters/year
ann_runoff = zeros(sim_len_yr ,1);   % meters/year
ann_overflow = zeros(sim_len_yr,1);   % meters/year
ann_del_water = zeros(sim_len_yr,1);   % meters/year

ann_del_C_del_t = zeros(sim_len_yr ,1); % ## change to ann_del_C_del_t

ann_WTD = zeros(sim_len_yr ,1);   % m  ## change to ann_wtd
ann_peat_water = zeros(sim_len_yr ,1);  % m3/m2  ## change to ann_peat_water
ann_total_water = zeros(sim_len_yr ,1);  % m3/m2  ## change to ann_total_water
ann_lagWTD = zeros(sim_len_yr ,1);  % years  ## change to ann_lag_wtd
ann_transmis = zeros(sim_len_yr ,1);   % relative hydraulic transmissivity (0-1)   
% annTHETA = zeros(sim_len_yr ,1);
% annWTD_VAR = zeros(sim_len_yr ,1);
ann_del_peatwater = zeros(sim_len_yr ,1); % ## change to ann_*
ann_net_water_in = zeros(sim_len_yr ,1);  % ## change to ann_*
age_depth = -9999 * ones(sim_len_yr ,30);
age_depth2 = -9999 * ones(sim_len_yr ,375);
growing_season_wtd = zeros(sim_len_yr,1);
num_growing_season_months = zeros(sim_len_yr,1);

% MONTHLY VECTORS BY TIME

mon_temp_forcing = zeros(sim_len_yr*12,1);    % degC
mon_precip_forcing = zeros(sim_len_yr*12,1);   % meters/month
mon_snowfall = zeros(sim_len_yr*12,1);   % meters/month
mon_rainfall = zeros(sim_len_yr*12,1);   % meters/month
mon_snowmelt = zeros(sim_len_yr*12,1);   % meters/month
mon_snowdepth = zeros(sim_len_yr*12,1);   % meters
mon_swe = zeros(sim_len_yr*12,1);   % meters
mon_snowsublimation = zeros(sim_len_yr*12,1);   % meters/month
mon_ALD = zeros(sim_len_yr*12,1);   % meters
mon_ALD2 = zeros(sim_len_yr*12,1);   % meters
mon_ALD3 = zeros(sim_len_yr*12,1);   % meters
mon_ALD1 = zeros(sim_len_yr*12,1);   % meters
mon_wtd = zeros(sim_len_yr*12,1);   % meters (positive down from surface)
mon_del_water = zeros(sim_len_yr*12,1);   % meters/month
mon_runon = zeros(sim_len_yr*12,1);   % meters/month
mon_runoff = zeros(sim_len_yr*12,1);   % meters/month
mon_transmis = zeros(sim_len_yr*12,1);   % dimensionless
mon_tot_water = zeros(sim_len_yr*12,1);   % meters
mon_overflow = zeros(sim_len_yr*12,1);   % meters/month
mon_aet = zeros(sim_len_yr*12,1);   % meters/month
mon_pet = zeros(sim_len_yr*12,1);   % meters/month
mon_water_content_error = zeros(sim_len_yr*12,1);   % meters

% for debugging:
mon_new_TOT_water_est = zeros(sim_len_yr*12,1);   % meters
mon_new_TOT_water_est2 = zeros(sim_len_yr*12,1);   % meters
mon_init_TOT_water_est = zeros(sim_len_yr*12,1);   % meters

mon_decompfact_water_top1000 = zeros(sim_len_yr*12,1);
mon_decompfact_temp_top1000 = zeros(sim_len_yr*12,1);

junk11 = zeros(sim_len_yr*3,2);
junk11_counter = 0;

ann_m_star_1yr = zeros(sim_len_yr, nveg); % Store mass remaining of each PFT after 1 year for comparison to litter bags
ann_m_star_2yr = zeros(sim_len_yr, nveg); % Store mass remaining of each PFT after 2 year for comparison to litter bags


% GIPL SOIL TEMP OUTPUT  % ** change to ann_snow_depth

ann_ALD1_max = zeros(sim_len_yr ,1);  % GIPL model annual max active layer depth to Tfr + FIT
ann_ALD2_max = zeros(sim_len_yr ,1);  % GIPL model annual max active layer depth to Tfr
ann_ALD3_max = zeros(sim_len_yr ,1);  % GIPL model annual max active layer depth to Tfr - FIT
ann_ALD_max = zeros(sim_len_yr ,1);  % selected one of above three to use
% soil_node_temp_month_save = zeros(sim_len_yr ,12,params_gipl.NumberOfSoilComputationNodes);
% soil_node_temp_month_saveMAT = zeros(sim_len_yr*12,params_gipl.NumberOfSoilComputationNodes);
ann_snow_depth = zeros(sim_len_yr ,1);  % GIPL model max annual snowdepth (meters?) %CT not actuually used.
ann_layer_Tmax = zeros(params_gipl.ndepth-1 ,2);   % stores annual max T of each soil layer for permmafrost of current year and previous year
tf_permafrost_layer_to10meters = zeros(params_gipl.ndepth-1,sim_len_yr);
zero_vec_SoilNodes = zeros(1,params_gipl.NumberOfSoilComputationNodes);

% vectors and arrays for the math

onevec = ones(sim_len_yr ,1);
epsvec = eps*ones(sim_len_yr ,1);
zerovec = zeros(sim_len_yr ,1);
onearr = ones(sim_len_yr ,nveg);
epsarr = eps*ones(sim_len_yr ,nveg);
topvec = zeros(1,nveg);
topval = 0;

% post-disturbance old and new carbon variables   % modified for HPM10PF

M_old = zeros(sim_len_yr ,1);        % = sum m across 'old' veg. types in cohort/layer i  
M_new = zeros(sim_len_yr ,1);        % = sum m across 'new' veg. types in cohort/layer i  
ann_M_old = zeros(sim_len_yr ,1);    % annual total 'old' M
ann_M_new = zeros(sim_len_yr ,1);    % annual total 'new' M
ann_NPP_old = zeros(sim_len_yr ,1);    % annual total 'old' NPP
ann_NPP_new = zeros(sim_len_yr ,1);    % annual total 'new' NPP
ann_resp_old = zeros(sim_len_yr ,1);    % annual total 'old' resp/decomp
ann_resp_new = zeros(sim_len_yr ,1);    % annual total 'new' resp/decomp

% ?? keep these??
cohort_age = zeros(sim_len_yr ,1);           % root-input-adjusted cohort age
cohort_age_temp = zeros(sim_len_yr ,1);           % root-input-adjusted cohort age
ann_resp_age = zeros(sim_len_yr ,1);   % age-mass weight of annual respiration
cohort_age2 = zeros(sim_len_yr ,nveg);           % root-input-adjusted cohort age
cohort_age_temp2 = zeros(sim_len_yr ,nveg);           % root-input-adjusted cohort age
ann_resp_age2 = zeros(sim_len_yr ,1);   % age-mass weight of annual respiration

% other useful things
start_day_of_month = [1,32,60,91,121,152,182,213,244,274,305,335];  % no leap years
end_day_of_month =  [31,59,90,120,151,181,212,243,273,304,334,365];  % no leap years
mid_day_of_month = round((start_day_of_month + end_day_of_month)/2);

% for testing
ann_snowsublimation_accum = 0;  % June 2018: Hamon PET method averages 60% of 10%_loss_per_month for Toolik temperatures 
WT_below_peat_counter = 0;

%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% READ IN MONTHLY WEATHER DRIVER and other input files
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% READ IN MONTHLY TEMPERATURE AND PRECIPITATION DRIVERS

%  header: years_forcing, months_forcing, precip_forcing_m, tsa_forcing_C 
%  years_forcing file (column #1) should be in years BP (1950 = 0 BP 1750 = 200 BP; 2010 = -60 BP)
%  months_forcing file (column #2) cycles through 1 to 12
%  precip_forcing_m (column #3) monthly precipitation in m/month liquid water
%  tsa_forcing_C file (column #4) surface air temperature in degrees C 
%  NOTE: climate forcing timespan should match inputs in parameter file, or should overwrite them
%  NOTE: climate forcing file should be chronological, so first line would have oldest driver (and largest year BP)

T = readtable(params.clim_in_name);
TA = table2array(T);

years_forcing = TA(:,1);
months_forcing = TA(:,2);
mon_precip_forcing = TA(:,3);
mon_temp_forcing = TA(:,4);

%  NOTE: climate forcing timespan should match inputs in parameter file, or should overwrite them
%  NOTE: climate forcing file should be chronological, so first line would have oldest driver (and largest year BP)
sim_start_check = max(years_forcing);  % = years_forcing(1);
sim_end_check = min(years_forcing);  % = years_forcing(end);

if (sim_start_check ~= params.sim_start)
    params.sim_start = params.sim_start 
    clim_file_start_year = sim_start_check
    stop
elseif (sim_end_check ~= params.sim_end)
    params.sim_end = params.sim_end 
    clim_file_end_year = sim_end_check
    stop
end

% accumulate monthly precip and average monthly temperature to annual values, just to have them
jan = 1:12:12*params.sim_len;
feb = jan + 1;
mar = jan + 2;
apr = jan + 3;
may = jan + 4;
jun = jan + 5;
jul = jan + 6;
aug = jan + 7;
sep = jan + 8;
oct = jan + 9;
nov = jan + 10;
dec = jan + 11;

ann_precip_forcing = mon_precip_forcing(jan) + mon_precip_forcing(feb) + ...
    mon_precip_forcing(mar) + mon_precip_forcing(apr) + mon_precip_forcing(may) + ...
    mon_precip_forcing(jun) + mon_precip_forcing(jul) + mon_precip_forcing(aug) + ...
    mon_precip_forcing(sep) + mon_precip_forcing(oct) + ...
    mon_precip_forcing(nov) + mon_precip_forcing(dec);

ann_temp_forcing = (mon_temp_forcing(jan) + mon_temp_forcing(feb) + ...
    mon_temp_forcing(mar) + mon_temp_forcing(apr) + mon_temp_forcing(may) + ...
    mon_temp_forcing(jun) + mon_temp_forcing(jul) + mon_temp_forcing(aug) + ...
    mon_temp_forcing(sep) + mon_temp_forcing(oct) + ...
    mon_temp_forcing(nov) + mon_temp_forcing(dec)) / 12.;


% READ IN ATMOSPHERIC DEL-14C AND INTERPOLATE TO ANNUAL

% read in annual atmos. Del-14C for 20000 BP to 2500 CE (-550 BP)
%  truncate to simulation years, and select RCP
%  del-14C data for RCPs to 2100 from Heather Graven (2015, PNAS, www.pnas.org/cgi/doi/10.1073/pnas.1504467112)
%      del-14C constant post-2100 (not realistic)

atm_del_c14_full_time_series = readtable(params.c14_in_name);  % Fix ?importdata? for UEF matlab??
atm_del_c14_years_input = atm_del_c14_full_time_series{:,1};  % curly brackets to extract values from table

% TRUNCATE DEL-C14 TO SIMULATION PERIOD
index_start = find(atm_del_c14_years_input == params.sim_start);
index_end = find(atm_del_c14_years_input == params.sim_end);
atm_del_c14_time_series = atm_del_c14_full_time_series(index_start:index_end,:);

ann_atm_del_c14_years = atm_del_c14_time_series{:,1};
ann_atm_del_c14 = atm_del_c14_time_series{:,1 + params.RCP_flag};
ann_atm_c14 = ann_atm_del_c14 / 1000. + 1;  

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% set up to save output for ?fancy graphs of M* and/or moss fraction of peat
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% *********** moss fraction *********************************************
% variables for binning M* and moss fraction of peat ********

flag_bins = 1;   %  If flag_bins = 1, compute M*, if 2, compute M* & moss_frac

if (flag_bins > 0.5) 

% ??? (next line)
    nbins = 250;                % for binning cohorts in output
% ??? (next line)
    maxheight = params.max_pot_peat_ht;   % total potential height (meters) 
    delx = maxheight/nbins;     % total possible ht (meters) ?? # of bins
    cohortheight = zeros(sim_len_yr,1);      % height of top of cohort above bottom of peat
    bin_M_star = 9999 * ones(nbins, sim_len_yr); % initialize array for M/M_0
end

if (flag_bins > 1.5) 
    mossfrac = zeros(sim_len_yr,1);          % initialize cohort mass fraction vector
    bin_moss_frac = -0.9999 * ones(nbins, sim_len_yr); % initialize array for fraction that is moss
end
 
% *****************************************************

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% INITIALIZE YEAR 1 COHORT, SOIL TEMPERATURE PROFILE, SWE, WFPS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% initialize surface cohort with aboveground litter inputs from all plant types

time(1) = 0.5 ;
time_month(1:12) = 0.5 + (0:(1/12):(1-1/12));

thick(1) = 0.05;  % placeholder value [meters] for first year NPP calculation

NPP = hpm20_mon_npp(params.wtd_0, params.wtd_0, params.ald_0, thick, params);  % use permafrost version

if (params.tf_old_new > 0.5)
    NPP = NPP  .* [old_new_ones old_new_zeros];   % modified for 2xN PFTs
end

m(1,:) = NPP .* params.ag_frac_npp;
m_0 = m;
m_star = m ./ (epsarr + m_0);
c14_m(1,:) = m(1,:) * ann_atm_c14(1);

M = sum(m,2);
M_0 = M;
M_star = M ./ (epsvec + M_0);
c14_M(1) = sum(m(1,:) .* c14_m(1,:)) / sum(m(1,:));

M_overlying(1) = 0;

ann_npp(1,:) = NPP;
ann_WTD(1) = params.wtd_0;
growing_season_wtd(1) = params.wtd_0;

% % Set initial soil temperatures
% T_profile_prev = params_gipl.T_init_gipl;
% % ??? (next line)  ?? add permafrost flag to parameters??
ann_ALD_max(1) = params.ald_0;  % arbitrary first year ALD (m)

% calculate layer density, thickness, and depth

dens = hpm20_mon_dens(M_star, M_overlying, params, onevec);

thick(1) = M(1) / (eps + dens(1));
zbottom(1) = thick(1);
prev_thick(1) = thick(1);
depth = cumsum(thick) - onevec * thick(1)/2;
ann_Z_total(1) = thick(1);
ann_M_total(1) = M(1);

tic;

flag1 = 0;   % set to 1 when simulation of dynamic water balance begins
% flag2 = 0;   % set to 1 when simulation of dynamic water balance begins
% ??? (next line)
profile_counter = 1;  % used to write out some profiles near end of run
accum_mon_del_water = 0;  % used to accumulate monthyly new water when it is too small to recompute WTD

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% LOOP THROUGH YEARS OF SIMULATION
%     NOTE: file ?hpm20_main_code_v1.m? has a number of climate perturbation simulation scripts
%           in the code just after the beginning of the annual loop (excluded here for now)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for iyear = 2: sim_len_yr 
    
    time(iyear) = iyear - 0.5;
    if (mod(iyear,500) == 0)   % tracks/writes out clock time per 500 y of simulation
        toc;
        tic;
        timex = (iyear - 0.5)
    end

%     COMPUTE DAILY T & P DRIVERS FOR GIPL2, USING 12 MONTHLY DRIVERS FOR THAT YEAR

    air_temp_month = mon_temp_forcing((iyear * 12  - 11) : (iyear * 12));
    precip_month = mon_precip_forcing((iyear * 12  - 11) : (iyear * 12));
    
    day_vec = 1:1:365;
%       day_vec = day_vec';

    daily_air_temp_for_year = interp1(mid_day_of_month, air_temp_month', day_vec, 'pchip', 'extrap');
    daily_precip_for_year = interp1(mid_day_of_month, precip_month, day_vec, 'pchip', 'extrap') / 30.5;
    % NOTE: daily precip will be daily drizzle/rain, without dry days, BUT ? 
    % NOTE: daily precip used in only GIPL2, and there only to generate snowpack, so drizzle is okay??

% zero out annual maximum values at beginning of year
    max_snow_depth = 0;  
    
    if (iyear == 2)   % NOTE: first year not simulated, so initialize things for start of year 2
        init_swe = 0.05;   % 5 cm snow water equivalent on 1 Jan of initial simulation year (year #2)
    end

    if (flag1 < 0.5)
        ann_WTD(iyear) = params.wtd_0;
        growing_season_wtd(iyear) = params.wtd_0;
%         mon_wtd_prev = params.wtd_0;

        Zstar = params.wfps_c1 * onevec + (params.wfps_c2 - params.wfps_c1)*((dens - params.min_bulk_dens)...
                     ./(dens - params.min_bulk_dens + params.wfps_c3));
        [wfps_init, PEAT_wat_init] = ...
                  hpm20_mon_wfps_only(params.wtd_0, Zstar,  ... 
                                             thick, depth, porosity, onevec, zerovec, params, params.ald_0);
        wfps = wfps_init;
    end
                                         
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%     LOOP THROUGH MONTHS OF YEAR
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    m_beginyear = m;

    for imonth = 1:1:12

        sim_month = (iyear-1) * 12 + imonth;
        time_month(sim_month) = 0.5 + (sim_month - 1)/12;
        
% extract month of daily T and P from annual (see above) for GIPL2 routine
        num_days_in_month = end_day_of_month(imonth) - start_day_of_month(imonth) + 1;
        daily_air_temp_for_month = daily_air_temp_for_year(start_day_of_month(imonth):end_day_of_month(imonth));
        daily_precip_for_month = daily_precip_for_year(start_day_of_month(imonth):end_day_of_month(imonth));

% COMPUTE MONTHLY SNOWFALL RAINFALL SNOWMELT, ETC.

        [month_snowfall, month_rainfall, month_snowmelt, month_snowdepth, month_snowsublimation, ...
            month_swe, ALFA, snowDepth, final_swe] = ...
            hpm20_mon_snowpack(num_days_in_month, daily_air_temp_for_month, daily_precip_for_month, init_swe, imonth, params_gipl);

% compute monthly snowsublimation in winter (if T_air < -1, snow sublimation = Hamon PET)
        [month_Hamon_pet] = hpm20_mon_pet(mon_temp_forcing(sim_month), imonth, num_days_in_month, params);
        
%        if (mon_temp_forcing(sim_month) < -1)
        if (final_swe > 0)
            month_snowsublim = min(final_swe, month_Hamon_pet);
            mon_snowsublimation(sim_month) = month_snowsublim;  % note, 19July2018, snow sublimation ~29% of annual snowfall for Lakkasuo
            final_swe = max(0, final_swe - month_snowsublim);
            ann_snowsublimation(iyear) = ann_snowsublimation(iyear) + month_snowsublim;
        end
            
        init_swe = final_swe;
        ann_rainfall(iyear) = ann_rainfall(iyear) + month_rainfall;
        ann_snowfall(iyear) = ann_snowfall(iyear) + month_snowfall;
        ann_snowmelt(iyear) = ann_snowmelt(iyear) + month_snowmelt;
%        ann_snowsublimation(iyear) = ann_snowsublimation(iyear) + month_snowsublimation;  # moved out of snowpack routine (see above)
        ann_max_snow_depth(iyear) = max(ann_max_snow_depth(iyear), month_snowdepth);

        mon_snowfall(sim_month) = month_snowfall;
        mon_snowmelt(sim_month) = month_snowmelt;
        mon_swe(sim_month) = month_swe;
        mon_rainfall(sim_month) = month_rainfall;
%        mon_snowsublimation(sim_month) = month_snowsublimation;
        mon_snowdepth(sim_month) = month_snowdepth;

        % for testing
                ann_snowsublimation_accum = ann_snowsublimation_accum + month_snowsublimation;
% CALL GIPL2 MONTHLY (save final day T(z) & SWE for next month; save T(z), SWE, rainfall, snowmelt, snowdepth)

%% call GIPL2 monthly 
%         [soil_layer_temp_month, soil_node_temp_month, soilTemp] = ...
%             hpm20_mon_gipl2(iyear, imonth, num_days_in_month, T_profile_prev, daily_air_temp_for_month, ...
%                     snowDepth, ALFA, depth, thick, wfps, porosity, dens, params, params_gipl);

%  calculate soil_layer_temperature for the month
    soil_node_temp_month = soil_node_temp_month_mat(sim_month, :);
    soil_layer_temp_month = 0.5*(soil_node_temp_month(1:params_gipl.ndepth-1) + soil_node_temp_month(2:params_gipl.ndepth));
% then go on with the calcs.
          layer_frozen = -((soil_layer_temp_month > 0) - 1);  % =0 if unfrozen, =1 otherwise (for transmissivity and infiltration and AET?)
        if (imonth == 1) 
            ann_layer_Tmax(:,2) = ann_layer_Tmax(:,1);
            ann_layer_Tmax(:,1) = soil_layer_temp_month;
        else
            ann_layer_Tmax(:,1) = max(ann_layer_Tmax(:,1), soil_layer_temp_month');
        end 
        
        if (imonth == 12)
            tf_permafrost_layer_to10meters(:,iyear) = (ann_layer_Tmax(:,1) < (params_gipl.Tfr + 0*params_gipl.FIT)) .* ...
                (ann_layer_Tmax(:,1) < (params_gipl.Tfr + 0*params_gipl.FIT));
        end
        
          
%% COMPUTE MONTHLY ALD  (invoke this only with permafrost, or use for seasonal frost layer as well?)
%   (invoke this only with permafrost during thaw season?). Montly ALD is used in water
%   balance, annual ALD is used for NPP.

        if (params.pf_flag > 0.5)   % params.pf_flag: 1=potential permafrost, 0=no permafrost
           if (iyear > 3)  % potential thaw months only
%            if (imonth > 3 && imonth < 11)  % potential thaw months only

%                [ALD1, ALD2, ALD3] = hpm20_mon_activelayer1(soil_node_temp_month, soil_layer_temp_month, params_gipl);
               [ALD1, ALD2, ALD3] = hpm20_mon_activelayer2(soil_node_temp_month_mat((sim_month -23):sim_month, 1:63), params_gipl);
                
%  which is best ALD metric, 1 or 2 or 3?
                mon_ALD1(sim_month) = ALD1;
                mon_ALD2(sim_month) = ALD2;
                mon_ALD3(sim_month) = ALD3;
                mon_ALD(sim_month) = min([ALD1, ALD2, ALD3]); %CT changed from mean
                
                ann_ALD1_max(iyear) = max(ALD1, ann_ALD1_max(iyear));
                ann_ALD2_max(iyear) = max(ALD2, ann_ALD2_max(iyear));
                ann_ALD3_max(iyear) = max(ALD3, ann_ALD3_max(iyear));

                ann_ALD_max(iyear) = max(mon_ALD(sim_month), ann_ALD_max(iyear));
           else
               mon_ALD(sim_month) = params.ald_0;
               ann_ALD_max(iyear) = params.ald_0;
            end   % % thaw months only
        end   % pf_flag = true (=1)


%% COMPUTE TERMS MONTHLY WATER BALANCE (delWater = Rain + snowmelt + runon ? AET ? runoff)

        if (flag1 > 0.5)  % compute monthly water balance
            
            [MonRunOff MonTransmis] = hpm20_mon_runoff(mon_wtd_prev, thick, depth, dens, onevec, layer_frozen, params,mon_ALD(sim_month));
            [MonRunOn] = hpm20_mon_runon(mon_wtd_prev, thick, depth, layer_frozen, params);
            [MonAET, MonPET] = hpm20_mon_aet(mon_wtd_prev, month_snowdepth, layer_frozen, mon_temp_forcing(sim_month), imonth, num_days_in_month, params);

%            mon_aet(sim_month) = MonAET;
            mon_aet(sim_month) = max( 0 , (MonAET - mon_snowsublimation(sim_month)));
            ann_aet(iyear) = ann_aet(iyear) + mon_aet(sim_month);
            mon_pet(sim_month) = MonPET;
            ann_pet(iyear) = ann_pet(iyear) + MonPET;
            mon_runoff(sim_month) = MonRunOff;
            mon_transmis(sim_month) = MonTransmis;
            ann_runoff(iyear) = ann_runoff(iyear) + MonRunOff;
            mon_runon(sim_month) = MonRunOn;
            ann_runon(iyear) = ann_runon(iyear) + MonRunOn;
            ann_transmis(iyear) = ann_transmis(iyear) + MonTransmis/12;

            mon_water_in = mon_rainfall(sim_month) + mon_snowmelt(sim_month) + mon_runon(sim_month);
            mon_water_out = mon_aet(sim_month) + mon_runoff(sim_month);
            mon_del_water(sim_month) = mon_water_in - mon_water_out;
            ann_del_water(iyear) = ann_del_water(iyear) + mon_water_in - mon_water_out;

            
% COMPUTE WTD AND WFPS(Z), only if monthly delta water > threshold

            accum_mon_del_water = accum_mon_del_water + mon_del_water(sim_month);

            if (abs(accum_mon_del_water) > params.del_water_threshold)  

                new_water = accum_mon_del_water;
                accum_mon_del_water = 0;
                mon_wtd_prev = mon_wtd(sim_month-1);

                Zstar = params.wfps_c1 * onevec + (params.wfps_c2 - params.wfps_c1)*((dens - params.min_bulk_dens)...
                     ./(dens - params.min_bulk_dens + params.wfps_c3));
                 
                WT_BELOW_PEAT_counter = WT_below_peat_counter;  % for tracking occurrence of deep water tables
            
                if (params.pf_flag > 0.5)   % (1=potential permafrost, 0=no permafrost)
                    [wfps, new_WTD, new_PEAT_water, wat_cont_error, new_TOT_water, WT_below_peat_counter, ...
                        new_TOT_water_est, new_TOT_water_est2] = ...
                          hpm20_mon_wtd_wfps(mon_wtd_prev, TOT_water, new_water, Zstar, dens, thick,   ... 
                                             depth, porosity, onevec, zerovec, params, mon_ALD(sim_month), WT_BELOW_PEAT_counter);
                else
                    ALD_no_pf = 100.;  % if no permafrost, set active layer to 100 m (i.e., ignore it)
                    [wfps, new_WTD, new_PEAT_water, wat_cont_error, new_TOT_water, WT_below_peat_counter, ...
                        new_TOT_water_est, new_TOT_water_est2] = ...
                          hpm20_mon_wtd_wfps(mon_wtd_prev, TOT_water, new_water, Zstar, dens, thick,   ... 
                                              depth, porosity, onevec, zerovec, params, ALD_no_pf, WT_BELOW_PEAT_counter);
                end
                
%                 function [new_wfps, new_WTD, new_PEAT_water, wat_cont_error, new_TOT_water, WT_below_peat_counter] = ...
%             hpm20_mon_wtd_wfps(WTD,init_TOT_water, new_water, Zstar, DENS,THICK,DEPTH,POROSITY, ...
%                                  ONEVEC,ZEROVEC,params ,ALT, WT_BELOW_PEAT_counter)

                mon_new_TOT_water_est(sim_month) = new_TOT_water_est;
                mon_new_TOT_water_est2(sim_month) = new_TOT_water_est2;
                mon_init_TOT_water_est(sim_month) = TOT_water + new_water;

                mon_wtd(sim_month) = new_WTD;
                mon_water_content_error(sim_month) = wat_cont_error;
                
                TOT_water = new_PEAT_water + (-new_WTD) * (new_WTD < 0);
                TOT_porosity = sum(thick .* porosity);
                ann_peat_water(iyear) = ann_peat_water(iyear) + new_PEAT_water/12;
                ann_total_water(iyear) = ann_total_water(iyear) + TOT_water/12;

                % DEBUG TEST                
%                 if (TOT_water < TOT_porosity && new_WTD < 0)
%                     junk = 1;
%                 end

% DEBUG TEST
             if (mon_del_water(sim_month) < 0)
                 junk = 1;
             end
            
             
            else  % very small net monthly water balance term -- skip computations
                mon_wtd(sim_month) = mon_wtd(sim_month-1);
                junk11_counter = junk11_counter + 1;
                junk11(junk11_counter,:) = [iyear imonth];

            end

% after computing new WTD, if it causes inundation > threshold, lose that excess water as runoff

            if (mon_wtd(sim_month) < -params.max_inundation)  %  params.max_inundation > 0; note: WTD < 0 if inundated.
                mon_overflow(sim_month) = -(mon_wtd(sim_month) - (-params.max_inundation));  % params.wtd_threshold
                ann_overflow(iyear) = ann_overflow(iyear) + mon_overflow(sim_month);
                mon_runoff(sim_month) = mon_runoff(sim_month) + mon_overflow(sim_month);
                ann_runoff(iyear) = ann_runoff(iyear) + mon_overflow(sim_month); 
                mon_wtd(sim_month) = -params.max_inundation;  
                TOT_water = TOT_water - mon_overflow(sim_month);
                ann_total_water(iyear) = ann_total_water(iyear) - mon_overflow(sim_month);
            end

            TOT_water_prev = TOT_water;

            mon_wtd_prev = mon_wtd(sim_month);
            mon_tot_water(sim_month) = TOT_water;

        else  % peat not deep enough to compute monthly water balance
            mon_wtd(sim_month) = params.wtd_0;
            wfps = wfps_init;
        end
 
        if (flag1 > 0.5)
            ann_WTD(iyear) = ann_WTD(iyear) + mon_wtd(sim_month)/12;
            growing_season_wtd(iyear) = growing_season_wtd(iyear) + mon_wtd(sim_month) * (mon_temp_forcing(sim_month) > 0);
            num_growing_season_months(iyear) = num_growing_season_months(iyear) + (mon_temp_forcing(sim_month) > 0);
        end
        
                
%%          COMPUTE MONTHLY DECOMPOSITION (?? Lose peat mass, but don?t recompute peat profile??)

        decompfact_water = hpm20_mon_decomp(depth, mon_wtd(sim_month), wfps, params, onevec, epsvec);
    
        decompfact_temp = (soil_layer_temp_month > (params_gipl.Tfr)) .* ... 
                                  (2 + 3*exp(-soil_layer_temp_month/9)).^((soil_layer_temp_month - 10)/10);
        cohort_decompfact_temp = interp1(params_gipl.soilNodeDepth(1:params_gipl.ndepth-1), decompfact_temp, depth);
                              
        mon_decompfact_water_top1000(sim_month) = mean(decompfact_water(1:min(params.sim_len, 1000)));
        mon_decompfact_temp_top1000(sim_month) = mean(cohort_decompfact_temp(1:min(params.sim_len, 1000)));
        k_mon = (decompfact_water .* cohort_decompfact_temp(:,1) * params.k_0 / 12) .* m_star;

        m_old = m;        
        m = m .* (onearr - k_mon);
        m_star = m ./ (m_0 + epsarr);
        c14_m = c14_m .* m ./ (m_old + epsarr);

        
%          SAVE MONTHLY RESULTS

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     END LOOP THROUGH MONTHS OF YEAR
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    end  % for loop imonth = 1:12

    if (flag1 > 0.5)
        growing_season_wtd(iyear) = growing_season_wtd(iyear) / num_growing_season_months(iyear);
    end

    ann_resp_C_array = (m_beginyear - m) / 2.;  % before new litter is added below
    
    %add litter tracker here for comparison to litter bags. Finds m/m0
    %after 12 and 24 months of decomposition
    if (iyear > 2)
        ann_m_star_1yr(iyear - 1, :) =  m(1, :) ./ (epsarr(iyear, :) + m_0(1,:)) + ...; %used calculation from later in code
                params.roots .* (m(2, :) ./ (epsarr(iyear, :) + m_0(2,:))); %1st year of roots go into 2nd litter cohort.
        ann_m_star_2yr(iyear -2, :) = (1- params.roots) .* m(2, :) ./ (epsarr(iyear, :) + m_0(2,:)) + ...;
                params.roots .* mean(m(3:4, :) ./ (epsarr(iyear, :) + m_0(3:4, :))); %because roots go in deeper: roots from 1,2,3 years old
    end
   
    
    
% sum up annual heterotrophic respiration by cohort and PFT; also 14C signal

    ann_RESP_C(iyear) = sum(sum(ann_resp_C_array));  % as carbon lost
    if (params.tf_old_new > 0.5)
        ann_RESP_C_new(iyear) = sum(sum(ann_resp_array(:,nveg/2+1:nveg),2));
        ann_RESP_C_old(iyear) = sum(sum(ann_resp_array(:,1:nveg/2),2));
    end

% move layers down one step to make room for new surface litter cohort (m,m_0,m_star,c14_m,depth)
% store masses in temporary vectors, adding zero to top, then removing final (bottom) zero to maintain same total size
% shift root input, and cohort masses m and m_0 and m_star down one cohort
             
    m_temp = [topvec; m];    % add zeros to top row
    m_0_temp = [topvec; m_0];
%    ms0agetemp = [topvec; m_0_age];
    m_star_temp = [topvec; m_star];
%    agebiastemp = [topval; age_bias];
    c14_mstemp = [topvec; c14_m];    % add zeros to top row
    depth_temp = [topval; depth];

    m_temp(end,:) = [];         % remove final row (of zeros) to maintain array size
    m_0_temp(end,:) = [];
%    ms0agetemp(end,:) = [];
    m_star_temp(end,:) = [];
%    agebiastemp(end) = [];
%    agebiastemparr = repmat(agebiastemp,1,nveg);
    c14_mstemp(end,:) = [];
    depth_temp(end) = [];
    depth = depth_temp;

% NEED TO FIGURE OUT THE NEXT LINE?S C14 COMPUTATION 
% (c14_mstemp = downshifted c14_m new cohort; mstemp is m at beginning of year; m_p12 is m at end of year)
    del_c14_ann_resp(iyear) = sum(sum(((c14_mstemp ./ (eps + m_temp) - 1) * 1000) .* (m_temp - m),2))/sum(sum((m_temp - m),2));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     COMPUTE ?ANNUAL? ALD AND WTD FOR NPP (?? Or instead compute monthly NPP within month loop??)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if iyear > params.lag_years     
    ann_lagWTD(iyear) = mean(ann_WTD(iyear:-1:(iyear-params.lag_years+1)));
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     COMPUTE ANNUAL NPP (or sum of monthly NPP)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% calculate annual productivity for each veg type as function of WTD

    NPP = hpm20_mon_npp(growing_season_wtd(iyear),ann_lagWTD(iyear), ann_ALD_max(iyear), thick, params);

 %    NPP = NPP * annTEMP_NPP_FACT(iyear);

    if (params.tf_old_new > 0.5) % if tf_old_new = 1, then implement old/new PFTs
           
        if (iyear < params.year_old2new)
            NPP = NPP .* [old_new_ones old_new_zeros];  % use 'old' PFTs 
        else
           NPP = NPP .* [old_new_zeros old_new_ones];  % use 'new' PFTs 
        end
        
    end

    ann_npp(iyear,:) = NPP(:);
    ann_NPP(iyear) = sum(NPP);
    if (params.tf_old_new > 0.5) 
        ann_NPP_old(iyear) = sum(NPP(1:nveg/2));
        ann_NPP_new(iyear) = sum(NPP(nveg/2+1:nveg));
    end

% determine root inputs
        
    rootin = hpm20_mon_rootin(depth, thick, params, NPP, growing_season_wtd(iyear), ann_ALD_max(iyear), ...?  
        ann_Z_total(iyear-1), onevec);

% for debugging?? compare these two root input values
    ann_ROOTIN(iyear) = sum(sum(rootin,2));
    ann_ROOTNPP(iyear) = sum(NPP .*params.bg_frac_npp);
    j5(iyear) = ann_ROOTIN(iyear) - sum(NPP .* params.bg_frac_npp); 

% move root inputs down one layer?
    rootin2 = [topvec; rootin];
    rootin2(end,:) = [];
    rootin = rootin2;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     ADD ABOVE-GROUND LITTER (NEW SURFACE COHORT) AND BELOW-GROUND LITTER
%         and decay 14C
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%  add roots

    m = m_temp + rootin;
    m_0 = m_0_temp + rootin;
    c14_m = c14_mstemp + rootin * ann_atm_c14(iyear);  % add root c14

%  add new top cohort
    
    m(1,:) = NPP .* params.ag_frac_npp;
    ann_AGMASSIN(iyear) = sum(ann_npp(iyear,:) .* params.ag_frac_npp);
    m_0(1,:) = m(1,:);
%    m_0_age(1,:) = m(1,:) * (params.sim_len - time(iyear));
%    age_bias(1) = 1;
    c14_m(1,:) = m(1,:) * ann_atm_c14(iyear);
    c14_m = c14_m * exp(-1/params.tau_c14);       % annual 14C loss to radioactive decay

    ann_resp_array = ann_resp_array + (m_temp - m) * 0.5 ;  % 0.5 factor for biomass to C
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%    UPDATE PEAT PROFILE (cohort m_star, M, m_0, m_star, bulk density, thickness)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    m_star = m ./ (epsarr + m_0);
    M = sum(m,2);
    M_tot = sum(M);
    prev_M_tot = M_tot;
    M_0 = sum(m_0,2);
%    M_0_age = sum(m_0_age,2);
    M_star = M ./ (epsvec + M_0);
    M_overlying = cumsum(M) - M;

    if (params.tf_old_new > 0.5) 
        M_old = sum(m(:,1:nveg/2),2);    % modified for old/new carbon
        ann_M_old(iyear) = sum(M_old);
        M_new = sum(m(:,nveg/2+1:nveg),2);    % modified for  old/new carbon
        ann_M_new(iyear) = sum(M_new);
    end
    
    del_M_tot(iyear) = M_tot - prev_M_tot; 

%  calculate cohort densities, thicknesses, and depths

    dens_temp = hpm_dens20(M_star, M_overlying, params, onevec);
    dens_prev_temp = [topval; dens_prev];
    dens_prev_temp(end) = [];
    dens_prev = dens_prev_temp;
    dens = max(dens_temp,dens_prev) .* (M > 0);   % don't let bulk density decline (??)
    dens_prev = dens;
    porosity = (onevec - dens/params.OM_dens) .* (M > 0);  % bulk density = mass peat ?? (volume organic matter + porosity)

%    dens_evolve(iyear,1) = dens(200);   %to see how densities of particular cohorts change with time
%    dens_evolve(iyear,2) = dens(500);
%    dens_evolve(iyear,3) = dens(1000);
%    dens_evolve(iyear,4) = dens(2000);
    
    prev_thick = thick;
    thick = M ./ (epsvec + dens);
    zbottom = cumsum(thick);
 
    if (max(zbottom) > params.start_depth && flag1 < 0.5)
        flag1 = 1
        dynamic_water_balance_year = iyear
        peat_height = max(zbottom)
        mon_wtd_prev = params.wtd_0
        TOT_water = sum(porosity .* thick .* wfps);
    end
    
    depth = (cumsum(thick) - thick/2) .* (M > 0);
%     if (iyear < params.sim_length)
%         depth(iyear+1:end) = 0;
%     end
    
    total_porosity = sum(thick .* porosity);   % total peat porosity = water content at saturation
    
    ann_Z_total(iyear) = depth(iyear)+thick(iyear)/2;    %    annZ_total(iyear) = sum(thick);   % = depth(iyear)+thick(1)/2;
    ann_del_Z_total(iyear) = ann_Z_total(iyear) - ann_Z_total(iyear-1);
    ann_M_total(iyear) = sum(M);

 % --------------------------
 % induce height-related shift to ombrotrophy (low O2 penetration/anoxia
 % scale length) + decrease in NPP
 %---------------------------
 if ann_Z_total(iyear) > x(2)
     params.anoxia_scale_length = x(5);
     params.NPP_rel = params.NPP_rel1 * x(3);
%  else
%      params.NPP_rel = params.NPP_rel1;
 end

 
    
% *********** moss fraction *********************************************
    % calculate moss fraction of peat in 'nbins' vertical bins over 'maxheight' meters from base 
    %       (can be greater than total peat height; missing value is
    %       -0.9999) (S. Frolking)
    if (flag_bins > 0.5)
        
        cohortheight = flipud(cumsum(flipud(thick))); 
        if (flag_bins > 1.5)
            mossfrac = (m .* (onevec * params.mosses)) ./ (M + eps);
        end

    
        x1 = 0.;
        for ix = 1:1:nbins
        
            if (x1 > max(cohortheight))
                break;
            end
            x2 = ix * delx;
    
            tf_bin = (cohortheight > x1) & (cohortheight <= x2);
            tf_bin_sum = sum(tf_bin);

            if(tf_bin_sum>0)
                bin_M_star(ix,iyear) = sum(M_star .* M .* tf_bin) / (sum(M .* tf_bin) + eps);
                if (flag_bins > 1.5)
                    bin_moss_frac(ix,iyear) = sum(mossfrac .* M .* tf_bin) / (sum(M .* tf_bin) + eps);
                end
            end
    
            x1 = x2;
        
        end     
          
    end
% *********************************************************************    

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%	?? RECOMPUTE WTD, OR NOT NECESSARY??
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%     SAVE ANNUAL RESULTS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END LOOP THROUGH YEARS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

end  % for loop of iyear = 2:sim_len_yr 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CALCULATE SOME FINAL METRICS AND WRITE OUT & PLOT RESULTS
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

age = time;
age_BP = time + + params.sim_end;
years_BP = flipud(time)+ params.sim_end;
M_TOTAL = sum(M);
M_TOTAL2 = sum(ann_del_C_del_t);
Z_TOTAL = depth(end);

k_mean = sum(m .* k_mon,2) ./ (M + epsvec);
c14_M = sum(c14_m,2) ./ (M + epsvec);
del_c14_M = (c14_M - 1) * 1000;   % final del-14C profile


wfps_c1a = 0.03;
wfps_c2a = 0.5;
wfps_c3a = 20;
zstar = params.wfps_c1 * onevec +(wfps_c2a - wfps_c1a)*((dens - params.min_bulk_dens)./(dens - (params.min_bulk_dens - wfps_c3a)));
sp_yld_profile = onevec - zstar + zstar .* ((onevec - zstar) / 0.01) .* exp(-max(0.5,depth)./zstar) .*  (onevec - exp(0.01*(onevec./zstar)));
sp_yld_profile =  max(zerovec,sp_yld_profile);   % specific yield profile

M_array = M * ones(1,nveg);
mfrac = m ./ M_array;
% 
% if (params.tf_old_new > 0.5)
%     annNPPmoss(:,1) = annNPP(:,1) + annNPP(:,6);  % modified for TOOLIK
% else
%     annNPPmoss(:,1) = annNPP(:,1);  % modified for TOOLIK
% end

annNPPmoss = sum(ann_npp .* (onevec * params.mosses), 2); 
annNPPvasc = sum(ann_npp,2) - annNPPmoss;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WRITE SUMMARY RESULTS TO SCREEN
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp(sprintf('total age (y): %d   total mass (kg C/m2): %d   total depth (m): %d',params.sim_len, M_TOTAL/2, Z_TOTAL));
disp(sprintf('total dC/dt (kg C/m2): %d ',M_TOTAL2));

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WRITE OUTPUT, AND SAVE WORKSPACE FOR PLOTTING
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% ***************************
% final core profile & processing

results_1 = [time depth M M_0 k_mean dens m del_c14_M];  

% data core ages, depths & find rows in results dataset
modCore_offset = min(params.obsCore_age) - params.sim_end;
modObs_ageIndex = params.obsCore_age - min(params.obsCore_age) + modCore_offset;
% checked against newly created vector Age_BP, which is derived from time
% vector that goes into results_1


rmseErr = rmse_modData(results_1(modObs_ageIndex, 2)*100, params.obsCore_depths.');
return
