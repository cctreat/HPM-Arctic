
% Revising GIPL2 soil thermal routine to shift from annual to monthly call;
% still operating on daily time step
%  S Frolking March 2017

% These lines above the function call will be moved out of the GIPL
% function into the main code


function [soil_layer_temp_month, soil_node_temp_month, swe_day_prev, soilTemp, max_snow_depth, ...
    mean_snow_depth, mean_swe, monthly_snowfall, monthly_rainfall, monthly_snowmelt] = ...
    hpm20_gipl2_daily(iyear, imonth, T_init_gipl, daily_air_temp_for_month, daily_precip_for_month, swe_day_prev, ...
    depth, thick, annWFPS, poros, dens, params, params_gipl)

% !!!!!!!! NOTE CHANGES TO FUNCTION CALL:
% ? CHANGE WFPS PROFILE INPUT IN FUNCTION CALL TO MONTHLY
% ? added imonth as input
% ? change ALD output to monthly
% ? keep max_snow_depth output?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
% uses UAF-GI GIPL2 model to compute daily soil temperature profiles, then averages to monthly
% GIPL2 code supplied by Sergei Marchenko (UAF), worked with version Claire Treat used in her MS
% 
num_days_per_month = 30;
day_vec = 1:1:num_days_per_month;
% day_vec = day_vec';

totalPeatDepth = max(depth);

daily_soil_temp_out = zeros(params_gipl.NumberOfSoilComputationNodes,num_days_per_month);
monthly_soil_temps = zeros(params_gipl.NumberOfSoilComputationNodes);
monthly_layer_temps = zeros(params_gipl.ndepth-1);
daily_soil_temps = zeros(params_gipl.NumberOfSoilComputationNodes,num_days_per_month);

% soil_temp_out = zeros(NumberOfSoilComputationNodes,num_months);

hpm_depth_index = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fpeat = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fwat = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fmin = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fporos = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fair = zeros(1,params_gipl.NumberOfSoilComputationNodes);
mean_ann_wfps = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cond_Th = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cond_Fr = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Th = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Fr = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Sl = zeros(1,params_gipl.NumberOfSoilComputationNodes);

U1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);
P1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Q1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);

T_air = zeros(1,num_days_per_month);
P = zeros(1,num_days_per_month);
snowfall = zeros(1,num_days_per_month);
rainfall = zeros(1,num_days_per_month);
snowsublimation = zeros(1,num_days_per_month);
snowmelt = zeros(1,num_days_per_month);
swe = zeros(1,num_days_per_month);
snowDepth = zeros(1,num_days_per_month);
ALFA = zeros(1,num_days_per_month);
zerovec = zeros(1,num_days_per_month);

    % see Dominik's paper for the unfrozen water function (negligible for peat, so maybe not a huge issue)

% **************************************
%  Inititialize peat and soil thermal properties each year
% **************************************

for (j=1:1:params_gipl.ndepth)

    if (totalPeatDepth >= params_gipl.soilLayerDepth(j))            % soil layer/node is above basal peat
        hpm_depth_index(j) = find(depth > params_gipl.soilLayerDepth(j), 1);
    else
        hpm_depth_index(j) = -99;
    end

    if (j == 1) 
        if (hpm_depth_index(j) > 0)  % peat
            
            mean_ann_wfps(j) = mean(annWFPS(1:hpm_depth_index(j)));
% !!!!!!!!CHANGE TO MONTHLY            
            mean_ann_wfps(j) = mean(annWFPS(1:hpm_depth_index(j)));
% !!!!!!!!
            fporos(j) = mean(poros(1:hpm_depth_index(j)));
            fpeat(j) = 1 - fporos(j);
            fmin(j) = 0;
            fwat(j) = mean(poros(1:hpm_depth_index(j))) .* annWFPS(1:hpm_depth_index(j));
% !!!!!!!!CHANGE TO MONTHLY            
            fwat(j) = mean(poros(1:hpm_depth_index(j))) .* WFPS_mon(1:hpm_depth_index(j));
% !!!!!!!!
            fair(j) = fporos(j) - fwat(j);
            
        else                         % sub-peat mineral soil
            
            fporos(j) = 0.65;     
            fpeat(j) = 0;
            fmin(j) = 1 - fporos(j);
            fwat(j) = 0.91 * fporos(j);    %  ORIGINALLY WAS 0.71 * fporos(j);
            fair(j) = fporos(j) - fwat(j);
            
        end

    else
        if (hpm_depth_index(j) > 0)  % peat
            
            mean_ann_wfps(j) = mean(annWFPS(hpm_depth_index(j-1):hpm_depth_index(j)));
% !!!!!!!!CHANGE TO MONTHLY            
            mean_mon_wfps(j) = mean(WFPS_mon(hpm_depth_index(j-1):hpm_depth_index(j)));
% !!!!!!!!
            fporos(j) = mean(poros(hpm_depth_index(j-1):hpm_depth_index(j)));
            fpeat(j) = 1 - fporos(j);
            fmin(j) = 0;
            fwat(j) = mean(poros(hpm_depth_index(j-1):hpm_depth_index(j)) .* ...
                       annWFPS(hpm_depth_index(j-1):hpm_depth_index(j)));
% !!!!!!!!CHANGE TO MONTHLY                
            fwat(j) = mean(poros(hpm_depth_index(j-1):hpm_depth_index(j)) .* ...
                       annWFPS(hpm_depth_index(j-1):hpm_depth_index(j)));
% !!!!!!!!
            fair(j) = fporos(j) - fwat(j);
            
        else                         % sub-peat mineral soil
            
            fporos(j) = 0.65;
            fpeat(j) = 0;
            fmin(j) = 1 - fporos(j);
            fwat(j) =  0.71 * fporos(j);
            fair(j) = fporos(j) - fwat(j);
            
        end
    end

    if (hpm_depth_index(j) > 0)  % peat
        
        Cond_Th(j) = 0.04 + 0.005 * mean_ann_wfps(j);   % from O'Donnell et al. 2009, via Wisser et al. 2011
        Cond_Fr(j) = 0.0141 + 0.0055 * mean_ann_wfps(j);
% !!!!!!!! CHANGE TO MONTHLY        
        Cond_Th(j) = 0.04 + 0.005 * mean_mon_wfps(j);   % from O'Donnell et al. 2009, via Wisser et al. 2011
        Cond_Fr(j) = 0.0141 + 0.0055 * mean_mon_wfps(j);
% !!!!!!!!
        Cvol_Th(j) = fpeat(j) * params_gipl.Cpeat + fwat(j) * params_gipl.Cwat;
        Cvol_Fr(j) = fpeat(j) * params_gipl.Cpeat + fwat(j) * params_gipl.Cice;  % UNFROZEN WATER CONTENT IN PEAT IS LOW (Wisser et al. 2011)
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);  % during phase change
        
    else                         % sub-peat mineral soil
        
        Cond_Th(j) = fmin(j) * params_gipl.Lmin + fwat(j) * params_gipl.Lwat + fair(j) * params_gipl.Lair;
        Cond_Fr(j) = fmin(j) * params_gipl.Lmin + fwat(j) * params_gipl.Lice + fair(j) * params_gipl.Lair;
        Cvol_Th(j) = fmin(j) * params_gipl.Cmin + fwat(j) * params_gipl.Cwat;
        Cvol_Fr(j) = fmin(j) * params_gipl.Cmin + (fwat(j) - params_gipl.Wunf) * params_gipl.Cice + params_gipl.Wunf * params_gipl.Cwat;  % UNFROZEN WATER CONTENT (= f(temp) in Wisser et al. 2011)
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);  % during phase change
        
    end

end

for j = params_gipl.ndepth+1:1:params_gipl.NumberOfSoilComputationNodes   % bedrock values
    
    if (params_gipl.soilLayerDepth(j) <= 30.)   % down to 30 meters
        
        Cvol_Fr(j) = 1.8e6;
        Cvol_Th(j) = 1.8e6;
        fwater(j) = 0.2;
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);
        Cond_Th(j) = 2.12;
        Cond_Fr(j) = 2.54;
%        fair(j) = 0;

    else                                           % below 30 meters
        
        Cvol_Fr(j) = 2.7e6;
        Cvol_Th(j) = 2.7e6;
        fwater(j) = 0.1;
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);
        Cond_Th(j) = 2.16;
        Cond_Fr(j) = 2.51;
%        fair(j) = 0;

    end
end

% if (iyear == 3550)   % for debugging?
%     pause;
% end

% !!!!!! MOVE SNOW CALCULATIONS OUT OF GIPL2
% ?????? MOVE SNOW CALCULATIONS OUT OF HPM?????
% **************************************
%  CALCULATED DAILY SNOW DEPTH, DENSITY, AND ALPHA VALUE
%     AND SAVE FOR DAILY TIME-STEP GIPL2 CALCULATIONS
% **************************************

max_snow_depth = 0;
monthly_snowfall = 0.;
monthly_rainfall = 0.;
monthly_snowmelt = 0.;
monthly_snowsublimation = 0;
mean_snow_depth = 0;
mean_swe = 0;

for iday = 1:1:num_days_per_month
     
    T_air(iday) = daily_air_temp_for_month(iday);  % °C
    P(iday) = daily_precip_for_month(iday);  % mm/mo to m/mo
    snowfall(iday) = P(iday) * (T_air(iday) < -1);%  based on Willmott et al. 1995
    rainfall(iday) = P(iday) - snowfall(iday);

%        T_air(imonth) = monthly_weather_time_series.data(iyear,imonth+1);  % °C
%        P(imonth) = monthly_weather_time_series.data(iyear,imonth+13)/1000;  % mm/mo to m/mo
%        snowfall(imonth) = P(imonth) * (T_air(imonth) < -1);%  based on Willmott et al. 1995
%        rainfall(imonth) = P(imonth) - snowfall(imonth); 

%        if (T_air_month < -1) %  based on Willmott et al. 1995
%            snowfall_month = P_month;
%           rainfall_month = 0;
%        else
%            snowfall_month = 0;
%            rainfall_month = P_month;
%        end

    snowsublimation(iday) = 0.1 * swe_day_prev / 30.5;  %  10 percent per month ???
    if (T_air(iday) > 1)
        snowmelt(iday) = (2.63 + 2.55*T_air(iday) + 0.0192*T_air(iday)*P(iday)) / 1000; % Willmott et al. 1995
    else
        snowmelt(iday) = 0;
    end
    swe(iday) = max(0, swe_day_prev + snowfall(iday) - snowmelt(iday) - snowsublimation(iday));
    swe_day_prev = swe(iday);
    if (swe(iday) > 0)
        snowDepth(iday) = swe(iday) * 1000 / params_gipl.snowDensity(imonth); % 1 m water = 1000 kg
        ALFA(iday) = 1.0/params_gipl.ALFA0 + snowDepth(iday)/(0.018 + 0.00087 * params_gipl.snowDensity(imonth));
    else
        snowDepth(iday) = 0;
        ALFA(iday) = params_gipl.ALFA0;
    end

    monthly_snowfall = monthly_snowfall + snowfall(iday);
    monthly_rainfall = monthly_rainfall + rainfall(iday);
    monthly_snowmelt = monthly_snowmelt + snowmelt(iday);
    monthly_snowsublimation = monthly_snowsublimation + snowsublimation(iday);
    max_snow_depth = max(max_snow_depth,snowDepth(iday));
    mean_snow_depth = mean_snow_depth + snowDepth(iday)/num_days_per_month;
    mean_swe = mean_swe + swe(iday)/num_days_per_month;

end  % loop through days of one month for snowpack calculations

% !!!!!! END SNOW CALCULATIONS


% call/implement GIPL2 model

% ********************************************
% LOOP THROUGH 1 MONTH OF GIPL2 SIMULATION
% ********************************************
 
for iday = 1:1:num_days_per_month
    
    maxABS = params_gipl.max_ABS;    % numerical iteration threshold

    soilTemp = zeros(1, params_gipl.NumberOfSoilComputationNodes);
    U1 = zeros(1, params_gipl.NumberOfSoilComputationNodes);
%        snowDepth(t)=0; edit 1-26-10
   
% 	within a loop --- time step = 1 (initial conditions)
    if  (iday == 1) % START WITH INITIAL CONDITIONS
            
        prevSoilTemp = T_init_gipl;
        U1 = prevSoilTemp;
   
    else   % not initial timestep

        U1 = prevSoilTemp;

% ---------------   
        S = 24*60*60;      %//   ! 1 day time step in seconds
%            S = 30.5*24*60*60; %//   ! 1 month (30.5 day) time step in seconds
        iter = 0; 
        
        while (iter < params_gipl.iter0 && maxABS > params_gipl.E0) %  (iter0 = 21, MAXIMUM?)
            
% 	---- computation of surface boundary coefficients  G1,G2
%            L0 = soilThermalConductivity(soilLayerDepth(1),U1(1));
%            L1 = soilThermalConductivity(soilLayerDepth(2),U1(2)); 	 
            if (U1(1) < (params_gipl.Tfr - params_gipl.FIT))
                L0 = Cond_Fr(1);
            elseif (U1(1) > (params_gipl.Tfr + params_gipl.FIT))
                L0 = Cond_Th(1);
            else
                L0 = 0.5 * (Cond_Th(1) + Cond_Fr(1));
            end
            if (U1(2) < (params_gipl.Tfr - params_gipl.FIT))
                L1 = Cond_Fr(2);
            elseif (U1(2) > (params_gipl.Tfr + params_gipl.FIT))
                L1 = Cond_Th(2);
            else
                L1 = 0.5 * (Cond_Th(2) + Cond_Fr(2));
            end

            H0 = params_gipl.soilLayerDepth(2) - params_gipl.soilLayerDepth(1);

            if snowDepth(iday) < params_gipl.E0 || prevSoilTemp(1) > params_gipl.E0 % no snow OR thawed surface%
%                if snowDepth(imonth) < E0 || prevSoilTemp(1) > E0 % no snow OR thawed surface
                G1 = 0.0;
                G2 = daily_air_temp_time_series_linear(iday);
%                G2 = monthly_weather_time_series.data(iyear,imonth+1);
                  
            elseif prevSoilTemp(1) <= 0.0 &&  snowDepth(imonth) < params_gipl.E0 % no snow and frozen surface
                G1 = 0.0;
                G2 = daily_air_temp_time_series_linear(iday);
%                G2 = monthly_weather_time_series.data(iyear,imonth+1);
 	                      
            else   % snowpack exists
    
%                ALFA = snowProperties(snowDensity, snowDepth(t));
                ALPHA = 1 / ALFA(iday);
%                ALPHA = 1 / ALFA(imonth);
%                C1 = heatCapacityDynWT(soilLayerDepth(1),prevSoilTemp(1));
                if (prevSoilTemp(1) < (params_gipl.Tfr - params_gipl.FIT))
                    C1 = Cvol_Fr(1);
                elseif (prevSoilTemp(1) > (params_gipl.Tfr + params_gipl.FIT))
                    C1 = Cvol_Th(1);
                else
                    C1 = Cvol_Sl(1);
                end
                W1 = 0.5 * (L0+L1);
                W2 = H0 * ALPHA / W1;
                W1 = 0.5 * power(H0,2) * C1 / W1 / S;
                G1 = 1.0 + W1 + W2;
                G2 = (W2 * daily_air_temp_time_series_linear(iday) + W1 * prevSoilTemp(1)) / G1;
%                G2 = (W2 * monthly_weather_time_series.data(iyear,imonth+1) + W1 * prevSoilTemp(1)) / G1;
                G1 = 1 / G1;
            end    

% 	----- Permutation and forward elimination
            P1(2) = G1;
            Q1(2) = G2;
                   
            for i = 2:1:(params_gipl.NumberOfSoilComputationNodes-1)
%                C1 = heatCapacityDynWT(soilLayerDepth(i),prevSoilTemp(i));   
                if (prevSoilTemp(i) < (params_gipl.Tfr - params_gipl.FIT))
                    C1 = Cvol_Fr(i);
                elseif (prevSoilTemp(i) > (params_gipl.Tfr + params_gipl.FIT))
                    C1 = Cvol_Th(i);
                else
                    C1 = Cvol_Sl(i);
                end
%                L2 = soilThermalConductivity(soilLayerDepth(i+1),prevSoilTemp(i+1));
                L2 = Cond_Fr(i+1) * (U1(i+1) < (params_gipl.Tfr - params_gipl.FIT)) + ... 
                      Cond_Th(i+1) * (U1(i+1) > (params_gipl.Tfr + params_gipl.FIT)) +  ...
                      0.5 * (Cond_Th(i+1) + Cond_Fr(i+1)) * ... 
                      ((U1(i+1) >= (params_gipl.Tfr - params_gipl.FIT)) && (U1(i+1) <= (params_gipl.Tfr + params_gipl.FIT)));

                H1 = params_gipl.soilLayerDepth(i+1) - params_gipl.soilLayerDepth(i);
                H2 = 0.5 * (H0 + H1);
                A1 = 0.5 * (L0 + L1) * S / C1 / (H0 * H2);
                B1 = 0.5 * (L1 + L2) * S / C1 / (H1 * H2);
                C0 = 1.0 + A1 + B1;
                P1(i+1) = B1 / (C0 - A1 * P1(i));
                Q1(i+1) = (A1 * Q1(i) + prevSoilTemp(i)) * P1(i+1) / B1;
                H0 = H1;
                L0 = L1;
                L1 = L2;
                      
%                heatCapacityOut(i,t) = C1;
%                thermalConductivityOut(i,t) = L2;
            end
%                      
% 	---- computation of the Lower boundary koef. G3 & G4
%            C1 = heatCapacityDynWT(soilLayerDepth(NumberOfSoilComputationNodes),...
%                     prevSoilTemp(NumberOfSoilComputationNodes));
            if (prevSoilTemp(params_gipl.NumberOfSoilComputationNodes) < (params_gipl.Tfr - params_gipl.FIT))
                C1 = Cvol_Fr(params_gipl.NumberOfSoilComputationNodes);
            elseif (prevSoilTemp(params_gipl.NumberOfSoilComputationNodes) > (params_gipl.Tfr + params_gipl.FIT))
                C1 = Cvol_Th(params_gipl.NumberOfSoilComputationNodes);
            else
                C1 = Cvol_Sl(params_gipl.NumberOfSoilComputationNodes);
            end
            G3 = 0.5 * power(H1,2) * C1 / L2 / S ;
            G4 = H1 * params_gipl.G0 + G3 * prevSoilTemp(params_gipl.NumberOfSoilComputationNodes);
            G3 = 1.0 / (1.0 + G3);
            G4 = G4 * G3;
                    
%  	---- Temperature computation in the last (deepest) grid node
            W1 = (G3 * Q1(params_gipl.NumberOfSoilComputationNodes) + G4) / ...
                   (1.0 - G3 * P1(params_gipl.NumberOfSoilComputationNodes));
% 	                  
            maxABS = abs(W1-U1(params_gipl.NumberOfSoilComputationNodes));
            U1(params_gipl.NumberOfSoilComputationNodes) = W1;
% 	                  
% 	 ---- Back substitution

            i = (params_gipl.NumberOfSoilComputationNodes-1);
            while (i>=1)
                W1 = P1(i+1) * U1(i+1) + Q1(i+1);
               
%	   ! check for the iterative convergence
                if (abs(W1-U1(i)) > maxABS)
                    maxABS = abs(W1-U1(i));
                end
                U1(i) = W1;
                           
                i = i - 1;

            end % END BACK-SUBSTITUTION WHILE (i>=1)                   
                       
            iter = iter + 1;
                    
        end % end while ((ITER < ITER0).AND.(maxABS > E0))
                            
        soilTemp = U1 ;

%         for (i=1:1:NumberOfSoilComputationNodes)%{//do i=1,N
        prevSoilTemp = soilTemp;                

        daily_soil_temp_out(1:params_gipl.NumberOfSoilComputationNodes,iday) = soilTemp;
%            soil_temp_out(1:NumberOfSoilComputationNodes, ((iyear - 1) * 12 + imonth)) = soilTemp';
%            soilT1 = soilTemp(1);
            
    end  % end of if loop for daily heat transfer computation
             

%  ---- Save daily temps for computing monthly values at end of year

    daily_soil_temps(:,iday) = soilTemp';
    
end  % loop of for loop through days (for iday = 1:num_days_per_month)

%  ---- soil layer temps equal average of the upper and lower nodal temps

daily_layer_temps = 0.5*(daily_soil_temps(1:params_gipl.ndepth-1,:) + daily_soil_temps(2:params_gipl.ndepth,:));
  
for jmonth=1:12   % 30 day months, starting with DOY 4
    monthly_layer_temps(:) = mean(daily_layer_temps(:,1:num_days_per_month,2));
    monthly_soil_temps(:) = mean(daily_soil_temps(:,1:num_days_per_month,2));
end

soil_layer_temp_month = monthly_layer_temps';  % IS THIS NEEDED?
soil_node_temp_month = monthly_soil_temps';   % IS THIS NEEDED?


% --- assemble soil temp output into array year, month, soil_level_above_bedrock

% monthly_soil_temp_array = zeros(num_years,12,ndepth-1);
% monthly_soil_temp_series = zeros(ndepth-1,num_years*12);
% for jyear = 1:1:num_years
%     for jmonth = 1:1:12
%         for jlevel = 1:1:ndepth-1
%             monthly_soil_temp_array(jyear,jmonth,jlevel) = mean(mean(soil_temp_out(jlevel:jlevel+1,((jyear-1)*366+(jmonth-1)*30+1):((jyear-1)*366+(jmonth)*30))));
%             monthly_soil_temp_series(jlevel,(jyear-1)*12+jmonth) = mean(mean(soil_temp_out(jlevel:jlevel+1,((jyear-1)*366+(jmonth-1)*30+1):((jyear-1)*366+(jmonth)*30))));
%         end
%     end
% end

        