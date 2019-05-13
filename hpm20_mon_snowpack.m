function [mon_snowfall, mon_rainfall, mon_snowmelt, mon_snowdepth, mon_snowsublimation, mon_swe, ALFA, snowDepth, final_swe] = ...
    hpm20_mon_snowpack(num_days_per_month,daily_air_temp_for_month, daily_precip_for_month, init_swe, imonth, params_gipl)

% compute daily snowdepth and ALFA for gipl2 model, and monthly snowpack stats
% rain/snow boundary and daily snowmelt from Willmott et al. 1995

snowfall = zeros(1,num_days_per_month);
rainfall = zeros(1,num_days_per_month);
snowsublimation = zeros(1,num_days_per_month);
snowmelt = zeros(1,num_days_per_month);
swe = zeros(1,num_days_per_month);
snowDepth = zeros(1,num_days_per_month);
ALFA = zeros(1,num_days_per_month);

mon_snowfall = 0;
mon_rainfall = 0;
mon_snowsublimation = 0;
mon_snowmelt = 0;
mon_swe = 0;
mon_snowdepth = 0;

swe_day_prev = init_swe;

for iday = 1:1:num_days_per_month
     
    snowfall(iday) = daily_precip_for_month(iday) * (daily_air_temp_for_month(iday) < -1);  
    rainfall(iday) = daily_precip_for_month(iday) - snowfall(iday);
%    snowsublimation(iday) = 0;  % move this out of function  # = 0.1 * swe_day_prev / 30.5;  %  10 percent per month; modify by daylength factor ???
    if (daily_air_temp_for_month(iday) > 1)
        snowmelt(iday) = max(2.63 + 2.55* daily_air_temp_for_month(iday) + ... 
            0.0192 * daily_air_temp_for_month(iday) * daily_precip_for_month(iday)) / 1000; 
        snowmelt(iday) = min(snowmelt(iday),swe_day_prev);
        
    else
        snowmelt(iday) = 0;
    end

    swe(iday) = max(0, swe_day_prev + snowfall(iday) - snowmelt(iday) - snowsublimation(iday));
    swe_day_prev = swe(iday);

    if (swe(iday) > 0)
        snowDepth(iday) = swe(iday) * 1000 / params_gipl.snowDensity(imonth); % 1 m water = 1000 kg/m2; snowdens range 200-500 kg/m3 % sensitivity analysis 2x snowpack
        ALFA(iday) = 1.0/params_gipl.ALFA0 + snowDepth(iday)/(0.018 + 0.00087 * ...            
                       params_gipl.snowDensity(imonth));
    else
        snowDepth(iday) = 0;
        ALFA(iday) = params_gipl.ALFA0;
    end

    mon_snowfall = mon_snowfall + snowfall(iday);
    mon_rainfall = mon_rainfall + rainfall(iday);
    mon_snowmelt = mon_snowmelt + snowmelt(iday);
    mon_snowsublimation = mon_snowsublimation + snowsublimation(iday);
    mon_snowdepth = mon_snowdepth + snowDepth(iday)/num_days_per_month; 
    mon_swe = mon_swe + swe(iday)/num_days_per_month;

end  % loop through days of one month for snowpack calculations

final_swe = swe(num_days_per_month);

