% script to read monthly soil node output for entire simulation from MAT file, and then generate
% some plots; needs 'soilLayerDepth' vector from parameter file

hpm20_mon_params_Seida;
% hpm20_mon_params_Lakkasuo;
params=load('hpm20_mon_param_vals');

nveg = params.num_veg;
sim_len_yr = params.sim_len

if (params.tf_old_new > 0.5)  % using Old-New carbon switch
    old_new_ones = ones(1,num_veg/2);
    old_new_zeros = zeros(1,num_veg/2);
end

if (params.pf_flag > 0.5) % simulation with permafrost
    hpm20_mon_gipl_params_pf;   % parameters for UAF GIPL 2.0 model with permafrost
else
    hpm20_mon_gipl_params_no_pf;   % parameters for UAF GIPL 2.0 model without permafrost
end

params_gipl = load('hpm20_mon_gipl_param_vals');   

soiltemp_filename = [params.out_name, '_soil_node_temp_month_save'];

load (soiltemp_filename);

% contains 'soil_node_temp_month_save': node temps in 3-d array (nyears, 12 months, nnodes)
% and contains 'SoilLayerDepth': node depths below surface (m) 

% extract layer temps (mean of adjacent nodes) for top 10 m (above bedrock)

array_dim = size(soil_node_temp_month_save);
nyears_tot = array_dim(1);
nnodes_tot = array_dim(3);

nnodes_above_bedrock = find(params_gipl.soilNodeDepth > 10,1) - 1;

% make array of monthly layer temps for final 'nyears' years of simulation

nyears = params.sim_len;
monthly_layer_temps_nyears = zeros(nyears*12,nnodes_above_bedrock-1);
for jyear = 1:1:nyears
    for jmonth = 1:1:12
         monthly_layer_temps_nyears((jyear-1)*12 + jmonth, :) = 0.5 * ...
            (soil_node_temp_month_save(nyears_tot - nyears + jyear,jmonth,1:nnodes_above_bedrock-1) + ...
             soil_node_temp_month_save(nyears_tot - nyears + jyear,jmonth,2:nnodes_above_bedrock));
    end
end

for j = 1:1:nyears/10
    monthly_layer_temps_nyears_10yr_mean(j,:) = sum(monthly_layer_temps_nyears(((j-1)*120 + 1):j*120,:))/120;
end

figure(41)
% contourf(monthly_layer_temps_nyears(1800:end,:)',15,':y')
% contourf(monthly_layer_temps_nyears(1:end,:)',15,':y')
%contourf(monthly_layer_temps_nyears_10yr_mean(1:end, 1:41)',15,':y') %top 3 m, decadal means
contourf(1:length(monthly_layer_temps_nyears_10yr_mean), soilNodeDepth(1:41), ...
    monthly_layer_temps_nyears_10yr_mean(1:end, 1:41)',15,':y') %top 3 m, decadal means
set(gca,'YDir','reverse') 

ystep = ( soilNodeDepth(1) + soilNodeDepth(2) ) / 2;
[X, Y] = meshgrid(1/12:1/12:nyears, ystep:ystep:10);

