function [new_wfps, new_WTD, new_PEAT_water, wat_cont_error, new_TOT_water, WT_below_peat_counter, ...
    new_TOT_water_est, new_TOT_water_est2] = ...
            hpm20_mon_wtd_wfps(WTD,init_TOT_water, new_water, Zstar, DENS,THICK,DEPTH,POROSITY, ...
                                 ONEVEC,ZEROVEC,params ,ALT, WT_BELOW_PEAT_counter)

% trying a new simplified version (May 2018) 
%  new_water equals monthly water change (m) (i.e., accum_mon_del_water)
%  guess WTD change as fraction of depth of new_water 
%  adjust estimate based on relative error of this first guess (i.e., interpolate or extrapolate?)
%  compute error (new_total_water vs TOT_water), make a final correction or add this to accum_mon_del_water for next month??

% OUTPUT VARIABLES
% ------------------
% new_wfps is updated peat corhort WFPS profile
% new_WTD is final water table depth below surface [m]
% new_peat_water is final water content of peat column (saturated + unsaturated) [m or m3/m2]
% wat_content_error is deviation of approximate solution to WTD for new water content [m or m3/m2]
% new_TOT_water is final (end of month) total water content (saturated + unsaturated + inundated) [m or m3/m2]
% initial (beginning of month) total water content (saturated + unsaturated + inundated) [m or m3/m2]
% WT_below_peat_counter is a counter for times when WT is below peat.

% INPUT VARIABLES
% ------------------
% WTD is initial water table depth below surface [m]
% init_TOT_water is initial (beginning of month) total water content (saturated + unsaturated + inundated) [m or m3/m2]
% new_water is monthly gain (>0) or loss (<0) of water (m or m3/m2)
% Zstar is factor for unsaturated WFPS, and is function of bulk density
% DENS is cohort bulk density profile [kg/m3]
% THICK is vector of cohort thicknesses (m)
% DEPTH is vector of cohort mid-point depths (m, positive down)
% POROSITY is vector of cohort porosities (- or m3/m3)
% ONEVEC & ZEROVEC are vectors of ones and zeros
% params is model input parameters
% ALT is current active layer thickness (is this needed?)
% WT_BELOW_PEAT_counter is a counter for times when WT is below peat.

% ------------------

% COMPUTED VARIABLES
% ------------------
% TOT_porosity is total peat column porosity  (- or m3/m3)
% zwtd is cohort distance above water table (used to compute WFPS of unsaturated peat) 
% wat_cont is cohort WFPS*porosity above the water table (m3/m3); below WT wat_cont = porosity
% min_peat_wat is peat profile column water if WT is at bottom of peat
% init_peat_water is peat profile column water at initial WTD (m or m3/m2)
% new_WTD_est is guess of new WTD position due to new_water addition (new_water can be negative)
% new_WFPS_est is WFPS profile associated with new_WTD_est
% new_TOT_water_est is peat profile water content associated with new_WFPS_est
% new_water_error is discrepancy between estimated water addition and new_water  
% new_WTD_est2 is correction to guess of new WTD position due to new_water addition (new_water can be negative)
% new_WFPS_est2 is WFPS profile associated with new_WTD_est2
% new_TOT_water_est2 is peat profile water content associated with new_WFPS_est2

% function parameter
wtd_step_fraction = 1.0;  % factor for change in WTD relative to depth of new_water

% ------------

TOT_porosity = sum(POROSITY .* THICK);
new_TOT_water = init_TOT_water + new_water;

%  first determine if WT is above peat surface, and set WTD and WFPS (WFPS=1)  and return
 
if (new_TOT_water >= TOT_porosity)  % peat is saturated, WTD at or above surface (i.e., WTD <= 0)
    
    new_WTD = -(new_TOT_water - TOT_porosity);
    new_wfps = ONEVEC;
    new_PEAT_water = TOT_porosity;
    wat_cont_error = 0;
    WT_below_peat_counter = WT_BELOW_PEAT_counter;
    new_TOT_water_est = new_TOT_water;
    new_TOT_water_est2 = new_TOT_water;

    return

end 

% ------------

%  second determine if total water would put WT below peat, and if so set WT in bottom peat layer, 
%    and compute WFPS and total water & return

zwtd = DEPTH - max(DEPTH);
zwtd = max(ZEROVEC, -zwtd);
wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
min_peat_water = sum((wat_cont .* THICK) .* POROSITY);

if (new_TOT_water < min_peat_water)
    new_WTD = max(DEPTH);
    new_wfps = wat_cont;
    new_PEAT_water = min_peat_water;
    wat_cont_error = -(new_TOT_water - new_PEAT_water);   % should be >0, meaning water error is gain
    new_TOT_water = new_PEAT_water;
    new_TOT_water_est = new_PEAT_water;
    new_TOT_water_est2 = new_PEAT_water;

    WT_below_peat_counter = WT_BELOW_PEAT_counter + 1;
    return;
end

% ------------

% Otherwise, take and guess at new WTD, and then approximate location based on error in guess  
%  NOTE that this could be put into a while loop that continued until the error was below some threshold

% first remove the inundation water, set WT to zero, and reduce water loss
% to that lost from peat itself.  NOTE: this condition should only be met
% if init_TOT_water > TOT_porosity AND(!) new_water < 0 AND(!) this loss of
% water (new_water<0) would bring new_TOT_water below TOT_porosity

if (init_TOT_water > TOT_porosity)
    new_water = new_water + (init_TOT_water - TOT_porosity);
    WTD = 0;
end
    
% estimate new WT position, keeping it at or below peat surface, and above the peat base
new_WTD_est = max(0, WTD - new_water * wtd_step_fraction);  
new_WTD_est = min(new_WTD_est, max(DEPTH));  

% determine how much new water the new WT position requires
zwtd = DEPTH - new_WTD_est;
zwtd = max(ZEROVEC, -zwtd);
new_WFPS_est = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
new_TOT_water_est = sum((new_WFPS_est .* THICK) .* POROSITY);
new_water_error = new_TOT_water_est - new_TOT_water;  

% adjust estimate of new WT position by the relative error in this new water estimate
% new_WTD_est2 = WTD - new_water * wtd_step_fraction * abs(new_water / (new_TOT_water_est + eps));
new_WTD_est2 = new_WTD_est + abs((new_WTD_est - WTD) / (new_water + eps)) *  new_water_error ;
new_WTD_est2 = min(new_WTD_est2,max(DEPTH));

% compute peat water content for this new WT position estimate
zwtd = DEPTH - new_WTD_est2;
zwtd = max(ZEROVEC, -zwtd);
new_WFPS_est2 = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
new_TOT_water_est2 = sum((new_WFPS_est2 .* THICK) .* POROSITY);
new_water_error2 = new_TOT_water_est2 - new_TOT_water;  

%  REPREAT CORRECTION A SECOND TIME
% adjust estimate of new WT position by the relative error in this new water estimate
% new_WTD_est2 = WTD - new_water * wtd_step_fraction * abs(new_water / (new_TOT_water_est + eps));
new_WTD_est3 = new_WTD_est2 + abs((new_WTD_est2 - WTD) / (new_water + eps)) * new_water_error2 ;
new_WTD_est3 = min(new_WTD_est3,max(DEPTH));

% compute peat water content for this new WT position estimate
zwtd = DEPTH - new_WTD_est3;
zwtd = max(ZEROVEC, -zwtd);
new_WFPS_est3 = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
new_TOT_water_est3 = sum((new_WFPS_est3 .* THICK) .* POROSITY);
% new_water_error3 = new_TOT_water_est3 - new_TOT_water;  

new_WTD = new_WTD_est3;
new_wfps = new_WFPS_est3;
new_PEAT_water = new_TOT_water_est3;
wat_cont_error = -(new_TOT_water - new_PEAT_water);   % make use of this error (e.g., pass along to next month)?
new_TOT_water = new_PEAT_water;

% SHOULD THERE BE NANOTHER REPEAT OF THE CORRECTION???

WT_below_peat_counter = WT_BELOW_PEAT_counter;  % do not increment

% for debugging (June 2018)
if (max(DEPTH) > 1) 
    junk = 1;
end
%  ADD any checks on WT position (e.g., is it within the peat?)

return;

% ------------



