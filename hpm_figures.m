%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PLOT FIGURES AND SAVE SOME TO FILES
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nyears = params.sim_len;

% ************************************
% FIGURE 1 - 4 panel final profiles
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(1)
subplot(1,4,1)
plot(m, depth,'LineWidth',2) 
set(gca,'YDir','reverse')
%ylabel('\fontsize{14}depth [m]')
xlabel('\fontsize{14}cohort mass [kg/m2]')
%title('\fontsize{14}cohort mass');
h1b=gca; 
set(h1b,'FontSize',14)

subplot(1,4,2)
plot(m_star, depth,M_star,depth,'LineWidth',2) 
set(gca,'YDir','reverse')
% ylabel('depth [m]')
title('\fontsize{14}cohort m-star values');
%             grs  minh wtms mins dshr hols lawn hums fthr ombs ombh evrs
%             trees
legend('\fontsize{10}PFT 1','\fontsize{10}PFT 2','\fontsize{10}PFT 3','\fontsize{10}PFT 4',...
    '\fontsize{10}PFT 5', '\fontsize{10}All',...
    'Location','SouthEast') 
% legend('\fontsize{10}grs','\fontsize{10}minh','\fontsize{10}mins','\fontsize{10}dshr',...
%     '\fontsize{10}wtms','\fontsize{10}hols','\fontsize{10}lawn','\fontsize{10}hums',...
%     '\fontsize{10}fthr','\fontsize{10}ombs','\fontsize{10}ombh','\fontsize{10}evrs',...
%     '\fontsize{10}trees','\fontsize{10}total','Location','SouthEast') 
h1c=gca; 
set(h1c,'FontSize',14)

subplot(1,4,3)
plot(dens, depth,'LineWidth',3) 
set(gca,'YDir','reverse')
% ylabel('depth [m]')
xlabel('\fontsize{14}bulk dens [kg/m3]');
%title('\fontsize{14}cohort m-star values');
h1d=gca; 
set(h1d,'FontSize',14)

subplot(1,4,4)
plot(hyd_trans_profile,depth, sp_yld_profile,depth, age/nyears,depth,'LineWidth',3)
set(gca,'YDir','reverse') 
title('\fontsize{14}final profile');
%xlabel('\fontsize{14}rel. trans.');
legend('\fontsize{14}rel.trans','\fontsize{14}sp.yield','\fontsize{14}age','Location','South')
h1e=gca; 
set(h1e,'FontSize',14)


h1c = gcf;
fig1name = [params.out_name, '_MstarPFT'];
saveas(h1c, fig1name,'jpg');

end

% *************************************
% FIGURE 2 - 4 panel summary time series
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(2)
subplot(4,1,1)
plot(time,ann_Z_total,age,-depth,'LineWidth',3)
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
legend('\fontsize{14}time-height','\fontsize{14}age-depth','Location','East')
ylabel('\fontsize{14}height or depth [m]')
xlabel('\fontsize{14}cohort or preatland age [y]')
% title('\fontsize{14}age-depthprofile #1');
h2a=gca; 
set(h2a,'FontSize',14)

subplot(4,1,2)
plot(time,ann_M_total,'LineWidth',3)
xlim([dynamic_water_balance_year+10 nyears+500]);
legend('\fontsize{14}time-mass','Location','East')
ylabel('\fontsize{14}total peat mass [kg/m2]')
xlabel('\fontsize{14}preatland age [y]')
% title('\fontsize{14}age-depthprofile #1');
h2b=gca; 
set(h2b,'FontSize',14)

subplot(4,1,3)
plot(time,ann_precip_forcing,time,-ann_WTD,'LineWidth',3)
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
% set(gca,'YDir','reverse')
legend('\fontsize{14}ann precip','\fontsize{14}WTD','Location','East')
ylabel('\fontsize{14}ann ppt or WTD [m]')
xlabel('\fontsize{14}peatland age [y]')
h2c=gca; 
set(h2c,'FontSize',14)
% title('\fontsize{14}age-depth profile #2');

subplot(4,1,4)
plot(years_BP,atm_del_c14_time_series{:,1},years_BP,del_c14_ann_resp,'LineWidth',3)
xlim([-100 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
set(gca,'XDir','reverse')
legend('\fontsize{14}atmosphere','\fontsize{14}respiration','Location','NorthWest')
ylabel('\fontsize{14}del-14C [o/oo]')
xlabel('\fontsize{14}years BP')
h2c=gca; 
set(h2c,'FontSize',14)
% title('\fontsize{14}age-depth profile #2');


end

% *************************************
% FIGURE 3 - 2 panel annual water flows
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(3)
subplot(2,1,1)
plot(ann_WTD,'LineWidth',3) 
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
set(gca,'YDir','reverse')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}WTD [m]')
h3a=gca; 
set(h3a,'FontSize',14)

subplot(2,1,2)
plot(time,ann_del_water, time,ann_precip_forcing, time,ann_aet, time,ann_runoff, time,ann_runon, time,-ann_WTD, time,-ann_ALD_max,'LineWidth',3) 
xlim([dynamic_water_balance_year+10 nyears+500]);
legend('\fontsize{14}dH2O_{cm}','\fontsize{14}ppt','\fontsize{14}ET','\fontsize{14}R_{off}','\fontsize{14}R_{on}','\fontsize{14}WT','\fontsize{14}ALT{cm}') 
legend('orientation','Horizontal','Location','North')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}depth [m]')
h3b=gca; 
set(h3b,'FontSize',14)

h3c = gcf;
fig3name = [params.out_name, '_watBal_diag'];
saveas(h3c, fig3name,'jpg');

end

% *************************************
% FIGURE 4 - 2 panel annual C flows
% *************************************

% plot or not

plot_flag = 0;

if (plot_flag > 1)

figure(4)
subplot(2,1,1)
plot(time, tot_npp/2, time, ann_RESP, time, del_C_del_t,'LineWidth',3)
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
legend('\fontsize{14}total NPP','\fontsize{14}ann resp','\fontsize{14}ann dC/dt') 
legend('orientation','Horizontal','Location','NorthEast')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}kgC/m2/y')
h4a=gca; 
set(h4a,'FontSize',14)

subplot(2,1,2)
plot(time, ann_NPPmoss,time,ann_NPPvasc,'LineWidth',2) 
%plot(time, annNPP,time,annROOTIN,'LineWidth',2) 
xlim([dynamic_water_balance_year+10 nyears+500]);
legend('\fontsize{10}moss','\fontsize{10}vascular')
%legend('\fontsize{10}grs','\fontsize{10}minh','\fontsize{10}mins','\fontsize{10}dshr',...
%    '\fontsize{10}wtms','\fontsize{10}hols','\fontsize{10}lawn','\fontsize{10}hums',...
%    '\fontsize{10}fthr','\fontsize{10}ombs','\fontsize{10}ombh','\fontsize{10}evrs',...
%    '\fontsize{10}total','Location','SouthEast') 
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}NPP [kg/m2/y]')
h4b=gca; 
set(h4b,'FontSize',14)

end

% *************************************
% FIGURE 5 - 3 panel final core profile by PFT
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(5)
subplot(1,3,1)
plot(M,time,'LineWidth',3)
set(gca,'YDir','reverse')
ylabel('\fontsize{14}cohort age [y]')
xlabel('\fontsize{14}cohort total mass [kg/m2]')
xlim([0 0.5])
h5a=gca; 
set(h5a,'FontSize',14)

subplot(1,3,2)
area(m)
view(90, 90)
legend('\fontsize{14}PFT 1','\fontsize{14}PFT 2','\fontsize{14}PFT 3','\fontsize{14}PFT 4',...
    '\fontsize{14}PFT 5',...
    'Location','SouthEast') 
%xlabel('\fontsize{14}cohort age [y]')
ylabel('\fontsize{14}cohort mass [kg/m2]')
ylim([0 0.25])
h5a=gca; 
set(h5a,'FontSize',14)

% subplot(2,1,2)
% area(rootin)
% ylabel('\fontsize{14}cohort m*')
% h6=gca; 
% set(h6,'FontSize',14)

subplot(1,3,3)
barh(depth,mfrac,'stack')
set(gca,'YDir','reverse')
ylabel('\fontsize{14}depth [m]')
xlim([0 1]);
xlabel('\fontsize{14}vegetation fraction of total mass');
% legend('\fontsize{10}grs','\fontsize{10}minh','\fontsize{10}mins','\fontsize{10}dshr',...
%    '\fontsize{10}wtms','\fontsize{10}hols','\fontsize{10}lawn','\fontsize{10}hums',...
%    '\fontsize{10}fthr','\fontsize{10}ombh','\fontsize{10}ombs','\fontsize{10}evrs',...
%    'Location','SouthEast') 
h5b=gca; 
set(h5b,'FontSize',14)

h5c = gcf;
fig5name = [params.out_name, '_pft_profile'];
saveas(h5c, fig5name,'jpg');

end

% *************************************
% FIGURE 6 - 2 panel C time series
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(6)
subplot(2,1,1)
plot(time,ann_Z_total, time, ann_M_total/100, age,-depth,'LineWidth',3)
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
legend('\fontsize{10}time-height','\fontsize{10}time-mass/100','\fontsize{10}age-depth','Location','East')
ylabel('\fontsize{14}[kg/m2] or [m]')
xlabel('\fontsize{14}cohort or peatland age [y]')
% title('\fontsize{14}age-depthprofile #1');
h6a=gca; 
set(h6a,'FontSize',14)

subplot(2,1,2)
plot(time,ann_NPP/2, time,annNPP_moss(:,1)/2, time, annNPP_herb/2, time, annNPP_woody/2,'LineWidth',2) 
% plot(time,ann_NPP/2, time,annNPP_moss(:,1)/2, time, annNPP_herb/2, time, annNPP_woody/2, time,ann_RESP_C, time,-ann_WTD,'LineWidth',2) 
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
legend('\fontsize{10}total NPP','\fontsize{10}moss NPP','\fontsize{10}vasc NPP','\fontsize{10}woody NPP')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}[kg C/m2/y]')
h6b=gca; 
set(h6b,'FontSize',14)

h6c = gcf;
fig6name = [params.out_name, '_H_M_NPP_R'];
saveas(h6c, fig6name,'jpg');

end

% *************************************
% FIGURE 7 - 2 panel annual water
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(7)
subplot(2,1,1)
plot(time,ann_del_water, time,ann_precip_forcing, time,ann_aet, time,ann_runoff, time,ann_runon, time,-ann_WTD, time,-ann_ALD_max,'LineWidth',3) 
xlim([dynamic_water_balance_year+10 nyears+500]);
legend('\fontsize{14}dH2O_{cm}','\fontsize{14}ppt','\fontsize{14}ET','\fontsize{14}R_{off}','\fontsize{14}R_{on}','\fontsize{14}WT','\fontsize{14}ALT{cm}') 
legend('orientation','Horizontal','Location','South')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}[m/y]')
h7b=gca; 
set(h7b,'FontSize',14)

subplot(2,1,2)
plot(ann_WTD,'LineWidth',3) 
xlim([dynamic_water_balance_year+10 nyears+500]);
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
set(gca,'YDir','reverse')
xlabel('\fontsize{14}time [y]')
ylabel('\fontsize{14}WTD [m]')
h7a=gca; 
set(h7a,'FontSize',14)

h7c = gcf;
fig7name = [params.out_name, '_water'];
saveas(h7c, fig7name,'jpg');

end

% *************************************
% FIGURE 8 - 2 panel core profile
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0)

figure(8)
% SF: removed plot of WTD and reconstructed water table depth (Aug. 2011)
% subplot(1,2,1)
  % semilogx(k_mean,depth,thick,depth,'LineWidth',3)
% plot(flipud(ann_WTD),depth, reconstWTD,depth,'lineWidth',3)
% set(gca,'YDir','reverse')
% hold on
% plot(zerovec,depth,'k','LineWidth',1)
% hold off
% set(gca,'YDir','reverse')
% legend('\fontsize{10}WTD','\fontsize{10}reconstr.','Location','SouthEast')
% ylabel('\fontsize{14}depth [m]')
% xlabel('\fontsize{14}WTD (positive down) [m]')
%plot(m_fast, depth) 
%set(gca,'YDir','reverse')
%ylabel('depth [m]')
%title('cohort fast-decomp mass');
% h8a=gca; 
% set(h8a,'FontSize',14)

subplot(1,2,1)
plot(dens/100,depth, hyd_trans_profile,depth, sp_yld_profile,depth, age/nyears,depth,'LineWidth',3)
set(gca,'YDir','reverse') 
title('\fontsize{14}final profile');
%xlabel('\fontsize{14}rel. trans.');
legend('\fontsize{10}bulkdens/100','\fontsize{10}rel.trans','\fontsize{10}sp.yield','\fontsize{10}age','Location','South')
h8b=gca; 
set(h8b,'FontSize',14)

subplot(1,2,2)
plot(del_c14_M,depth, flipud(atm_del_c14_time_series{:,1}),depth,'LineWidth',3)
set(gca,'YDir','reverse') 
title('\fontsize{14}final profile');
xlabel('\fontsize{14}del-14C [o/oo]')
legend('\fontsize{12}peat','\fontsize{12}atmosphere','Location','SouthEast')
h8c=gca; 
set(h8c,'FontSize',14)

h8c = gcf;
fig8name = [params.out_name, '_final_profile'];
saveas(h8c, fig8name,'jpg');
end

% ******************************************************************
% FIGURE 24 - moss fraction
% *************************************

% plot or not

plot_flag = -1;

if (plot_flag > 0 && flag_bins > 1.5)

[XXX,YYY] = meshgrid(1:1:itime, 1:1:nbins);
xxx = [1:1:itime];
yyy = [0,delx,(maxheight-delx)];
clims = [-0.1 1];

%tf_sp = bin_moss_frac == -0.9999;
% bin_moss_frac = bin_moss_frac.*(bin_moss_frac>=0)+bin_moss_frac.*(bin_moss_frac<0).*NaN;

figure(24)
% contourf(XXX, YYY, bin_moss_frac,100,'LineStyle','none')
% xlim([0 itime])
%imagesc(xxx,fliplr(yyy),bin_moss_frac,clims);
h=imagesc(xxx,yyy,bin_moss_frac,clims);
colormap([[1,1,1];jet]);

%set(h,'alphadata',~isnan(bin_moss_frac))

hold on;

plot(xxx,(ann_Z_total-ann_WTD),'linewidth',1.5,'Color',[0 0 0]);

set(gca,'ydir','normal');

title('\fontsize{30}Moss fraction of peat');
colorbar('off'); colorbar('location','eastoutside');
% ylim([0 nbins])
% zlim([0 1])
% colorbar
% caxis([0 1])
xlabel('\fontsize{20}Simulation time [yr]')
ylabel('\fontsize{20}Peat height [m]')
zlabel('\fontsize{20}Moss fraction of peat')
h24a=gca;
set(h24a,'FontSize',18)


h24c = gcf;
fig24name = [params.out_name, '_mossfrac'];
saveas(h24c, fig24name,'jpg');

end

% ******************************************************************
% FIGURE 25 - M_Star
% *************************************

% plot or not

plot_flag = 1;

if (plot_flag > 0 && flag_bins > 0.5)

[XXX,YYY] = meshgrid(1:1:iyear, 1:1:nbins);
% xxx = [1:1:itime];
xxx = [iyear:-1:1];
yyy = [0,delx,(maxheight-delx)];
%clims = [-0.1 1.7];
clims = [-0.1 1.0];

%tf_sp = bin_moss_frac == -0.9999;
% bin_moss_frac = bin_moss_frac.*(bin_moss_frac>=0)+bin_moss_frac.*(bin_moss_frac<0).*NaN;

figure(25)
% contourf(XXX, YYY, bin_moss_frac,100,'LineStyle','none')
% xlim([0 itime])
%imagesc(xxx,fliplr(yyy),bin_moss_frac,clims);
% log_bin_M_star = (log(bin_M_star) - log(min(min(bin_M_star))) + 0.001) .* (bin_M_star<3) + ...
%                    (bin_M_star>=3).*-0.999;

log_bin_M_star = (log10(bin_M_star) - log10(min(0.02,min(min(bin_M_star))))) .* (bin_M_star<3) + ...
                  (bin_M_star>=3) * -0.999;
%bin_M_star_plot = bin_M_star .* (bin_M_star<3) - 0.9999 * (bin_M_star>=3);

h=imagesc(xxx,yyy,log_bin_M_star,clims);
%h=imagesc(xxx,yyy,bin_M_star_plot,clims);
colormap([[1,1,1];jet]);

%set(h,'alphadata',~isnan(bin_moss_frac))

hold on;

plot(xxx,fliplr(ann_Z_total-ann_WTD),'linewidth',1.5,'Color',[0 0 0]); % [0 0 0] - black
if params.pf_flag > 0.5
    plot(xxx,(ann_Z_total-ann_ALD_max),'linewidth',0.5,'Color',[0 0 1]);  % [0 0 1] - blue
end

set(gca,'ydir','normal');

title('\fontsize{20}fraction mass remaining of peat');
colorbar('off'); cb = colorbar('location','eastoutside'); shading flat;

set(cb,'fontsize',12);
caxis_ticks = 10.^(get(cb,'YTick')) * min(0.02,min(min(bin_M_star)));
caxis_ticks_str = Num2CellStr(caxis_ticks);
% for jax = 1:1:6
%     if (caxis_ticks(jax) < 0.1)
%         cax(jax) = num2str(caxis_ticks(jax),3);
%     else
%         cax(jax) = num2str(caxis_ticks(jax),4);
%     end
% end
% caxis_all = cell(cax);

set(cb,'YTickLabel',caxis_ticks_str);
   
ylim([0,2]) %max(depth)])
% zlim([0 1])
% colorbar
% caxis([0 1])
xlabel('\fontsize{20}years BP')
ylabel('\fontsize{20}Peat height [m]')
% zlabel('\fontsize{20}fraction mass remaining of peat')
% h25a=gca;
% set(h25a,'FontSize',18)


h25c = gcf;
fig25name = [params.out_name, '_log_m_star'];
saveas(h25c, fig25name,'jpg');

% ******************************************************************

end

% *************************************
% FIGURE 9 - 1 panel CAR time series
% *************************************

% plot or not

plot_flag = 0;

if (plot_flag > 0)

figure(9)
% subplot(2,1,1)
plot(time,M*500, MB930_age,MB930_CAR,'LineWidth',2)
xlim([dynamic_water_balance_year+10 nyears+500]);
legend('\fontsize{12}HPM10','\fontsize{12}MB930','Location','North')
ylabel('\fontsize{14}[gC/m2/y]')
xlabel('\fontsize{14}cohort age [years BP]')
title('\fontsize{24}apparent C accumulation rate');
h9a=gca; 
set(h9a,'FontSize',14)

h9c = gcf;
fig9name = [params.out_name, '_CAR'];
saveas(h9c, fig9name,'jpg');

end

% *************************************
% FIGURE 10 - 1 panel cohort depth time series
% *************************************

% plot or not

plot_flag = -1;

if (plot_flag > 0)

figure(10)

plot(age_depth,'.','MarkerSize',3)
hold on
plot(time,ann_Z_total,'LineWidth',4,'Color','black')
hold off
% xlim([0 max(time)+500])
xlim([0 max(time)])                   % axis limit for HPM CAR paper
% ylim([0.001 max(ann_Z_total)+0.5])
ylim([0 5.5])                         % axis limit for HPM CAR paper       
ylabel('\fontsize{14}peat height [m]')
xlabel('\fontsize{14}peatland age [y]')
title('\fontsize{14}500-year interval peat cohort trajectories');
set(gca,'XDir','reverse') 
h10a=gca;  
set(h10a,'FontSize',14)

h10b = gcf;
fig10name = [params.out_name, '_cohort_depth'];
saveas(h10b, fig10name,'jpg');

figure(30)

plot(age_depth2,'.','MarkerSize',3)
hold on
plot(time,ann_Z_total,'LineWidth',3,'Color','black')
hold off
% xlim([0 max(time)+500])
xlim([8300 max(time)])                   % axis limit for HPM CAR paper
% ylim([0.001 max(ann_Z_total)+0.5])
ylim([4. 5.3])                         % axis limit for HPM CAR paper       
ylabel('\fontsize{14}peat height [m]')
xlabel('\fontsize{14}peatland age [y]')
title('\fontsize{14}25-year interval peat cohort trajectories');
set(gca,'XDir','reverse') 
h30a=gca;  
set(h30a,'FontSize',14)

h30b = gcf;
fig30name = [params.out_name, '_cohort_depth2a'];
saveas(h30b, fig30name,'jpg');

end

% *************************************
% FIGURE 11 - 1 panel age-depth profile
% *************************************

% ** +/- one sigma on MB930 core radiocarbon dates (x is age range in cal.yrsBP, y is depth in m)
MB930_x1 = [540 730]; MB930_y1 = [0.395 0.395];
MB930_x2 = [2782 3472]; MB930_y2 = [1.8 1.8];
MB930_x3 = [3469 3550]; MB930_y3 = [1.905 1.905];  % this one looks odd, 
MB930_x4 = [6313 6724]; MB930_y4 = [3.005 3.005];
MB930_x5 = [7792 9423]; MB930_y5 = [4.805 4.805];

% plot or not

plot_flag = 0;

if (plot_flag > 0)

figure(11)
% subplot(2,1,2)
plot(age,depth,'LineWidth',3)
set(gca,'YDir','reverse') 
hold on
plot(MB930_age_depth_time_series.data(:,2),MB930_age_depth_time_series.data(:,1),'--ko',...
     'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',11)
% plot(MB930_x1,MB930_y1,MB930_x2,MB930_y2, MB930_x3,MB930_y3, MB930_x4,MB930_y4,MB930_x5,MB930_y5,'-k','LineWidth',1) 
hold off
set(gca,'YDir','reverse') 
% ylim([0.001 max(ann_Z_total)+0.5])
legend('\fontsize{14}HPM','\fontsize{14}MB930','Location','NorthEast')
ylabel('\fontsize{14}peat depth [m]')
xlabel('\fontsize{14}peat age [cal. y BP]')
title('\fontsize{14}age-depth');
h11a=gca; 
set(h11a,'FontSize',14)

h11c = gcf;
fig11name = [params.out_name, '_age_depth'];
saveas(h11c, fig11name,'jpg');

end

% *************************************
% FIGURE 12 - 1 NEE vs. CAR time series
% *************************************

% generate 100-year averages of NEE and CAR
navg = floor(nyears / 100.);
NEE_100yr = zeros(navg,1);
CAR_100yr = zeros(navg,1);

for iavg = 1:1:navg
    NEE_100yr(iavg) = mean(ann_del_C_del_t(1+(iavg-1)*100:iavg*100));
    CAR_100yr(iavg) = mean(M(1+(iavg-1)*100:iavg*100)/2);  % divide by 2 for mass to C
end
NEE_100yr = flipud(NEE_100yr);
time_avg = 100*[1:1:navg] - 50;

% save this into a file
results_12 = [time_avg' NEE_100yr CAR_100yr];
fname12 = [params.out_name, '_100y_NCB_CAR.txt'];
fid12 = fopen(fname12,'w');  % time series of carbon dynamics

fprintf(fid12,'HPM10 output - 100-yr mean time series of carbon dynamics - units: NCB, CAR: kg C/m2/y; time: y \n');
fprintf(fid12,' time     NCB   CAR\n');
fprintf(fid12,'%7.1f %9.5f %9.5f \n', results_12');
status = fclose(fid12);

% plot or not

plot_flag = -1;

if (plot_flag > 0)

figure(12)
% subplot(2,1,2)
plot(time_avg,NEE_100yr,time_avg,CAR_100yr,'LineWidth',3)
hold on
plot(time,zerovec,'k','LineWidth',1)
hold off
% set(gca,'YDir','reverse') 
% ylim([0.001 max(ann_Z_total)+0.5])
legend('\fontsize{14}HPM12 NCB','\fontsize{14}HPM12 CAR','Location','North')
ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
xlabel('\fontsize{14}years BP')
title('\fontsize{14}simulated NCB and CAR: 100-yr averages');
h12a=gca; 
set(h12a,'FontSize',14)

h12c = gcf;
fig12name = [params.out_name, '_NEE_vs_CAR'];
saveas(h12c, fig12name,'jpg');

end

% *************************************
% FIGURE 14 - 19 - smoothed summary plots (14-16 full simulation; 17-19 from 1801 to end)
% *************************************

yr_1800_CE = 150 - min(years_BP);

% results_13_smooth: time  annTEMP   annPPT    annET   annRUNOFF  annWTD   annNPP    annDECOMP  annNCB    annNPP_herb annNPP_woody annNPP_moss peat_mass peat_height \n');

% plot or not

plot_flag = 1;

if (plot_flag > 0)

    figure(14)
    subplot(4,2,1)  % annual temperature (degC)
    plot(results_13_smooth(:,1),results_13_smooth(:,2),'LineWidth',3)
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    title('\fontsize{10}annual temperature (degC): 100-yr loess smoothing');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,2)   % annual precipitation (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,3),'LineWidth',3)
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    title('\fontsize{10}annual precipitation (m/y): 100-yr loess smoothing');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,3)  % annual ET (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,4),'LineWidth',3)
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    title('\fontsize{10}annual ET (m/y): 100-yr loess smoothing');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,4)  % annual runoff (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,5),'LineWidth',3)
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    title('\fontsize{10}annual runoff (m/y): 100-yr loess smoothing');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,5)  % annual WTD (m below surface)
    plot(results_13_smooth(:,1),results_13_smooth(:,6),'LineWidth',3)
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    set(gca,'YDir','reverse') 
    title('\fontsize{10}annual WTD (m): 100-yr loess smoothing');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,6)  % annual NPP and decomposition (kg C/m2/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,7),results_13_smooth(:,1),-results_13_smooth(:,8),'LineWidth',3)
    hold on
    plot(time,zerovec,'k','LineWidth',1)
    hold off
    % legend('\fontsize{10}NPP','\fontsize{10}decomposition','Location','NorthEast')
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{10}years BP')
    title('\fontsize{10}ann. NPP & decomp. (kg C/m2/y): 100-yr loess sm.');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,7)  % annual precipitation (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,9),'LineWidth',3)
    hold on
    plot(time,zerovec,'k','LineWidth',1)
    hold off
    % ylabel('\fontsize{10}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}ann. net C bal. (kg C/m2/y): 100-yr loess sm.');
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    subplot(4,2,8)  % annual peat height (m) and mass/100 (kg/m2)
    plot(results_13_smooth(:,1),0.02*results_13_smooth(:,13),results_13_smooth(:,1),results_13_smooth(:,14),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}ann. peat H (m) & M/100 (kg/m2): 100-yr loess sm.');
    legend('\fontsize{10}mass/100','\fontsize{10}height','Location','NorthWest')
    xlim([0 max(time)+100])
    % xlim([max(time)-400 max(time)])                   

    h14a=gca; 
    set(h14a,'FontSize',8)

    h14c = gcf;
    fig14name = [params.out_name, '_smoothed_outputs'];
    saveas(h14c, fig14name,'jpg');




    figure(18)
    subplot(2,2,1)  % annual temperature (degC)
    plot(results_13_smooth(:,1),results_13_smooth(:,2),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{14}years BP')
    title('\fontsize{10}annual temperature (degC): 100-yr loess smoothing');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,2)   % annual precipitation (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,3),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{14}years BP')
    title('\fontsize{10}annual precipitation (m/y): 100-yr loess smoothing');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,3)  % annual ET (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,4),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}annual ET (m/y): 100-yr loess smoothing');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,4)  % annual runoff (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,5),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}annual runoff (m/y): 100-yr loess smoothing');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    h18a=gca; 
    set(h18a,'FontSize',10)

    h18c = gcf;
    fig18name = [params.out_name, '_smoothed_P_T_ET_RO_final'];
    saveas(h18c, fig18name,'jpg');

    figure(19)
    subplot(2,2,1)  % annual WTD (m below surface)
    plot(results_13_smooth(:,1),results_13_smooth(:,6),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{14}simulation years (end = 2100 CE)')
    set(gca,'YDir','reverse') 
    title('\fontsize{10}annual WTD (m): 100-yr loess smoothing');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,2)  % annual NPP and decomposition (kg C/m2/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,7),results_13_smooth(:,1),-results_13_smooth(:,8),'LineWidth',3)
    hold on
    plot(time,zerovec,'k','LineWidth',1)
    hold off
    legend('\fontsize{10}NPP','\fontsize{10}decomposition','Location','SouthWest')
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    % xlabel('\fontsize{14}years BP')
    title('\fontsize{10}ann. NPP & decomp. (kg C/m2/y): 100-yr loess sm.');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,3)  % annual precipitation (m/y)
    plot(results_13_smooth(:,1),results_13_smooth(:,9),'LineWidth',3)
    hold on
    plot(time,zerovec,'k','LineWidth',1)
    hold off
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}ann. net C bal. (kg C/m2/y): 100-yr loess sm.');
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    subplot(2,2,4)  % annual peat height (m) and mass/100 (kg/m2)
    plot(results_13_smooth(:,1),0.02*results_13_smooth(:,13),results_13_smooth(:,1),results_13_smooth(:,14),'LineWidth',3)
    % ylabel('\fontsize{14}NCB or CAR [g C/m2/yr]')
    xlabel('\fontsize{10}simulation years (end = 2100 CE)')
    title('\fontsize{10}ann. peat H (m) & M/100 (kg/m2): 100-yr loess sm.');
    legend('\fontsize{10}mass/100','\fontsize{10}height','Location','South')
    % xlim([0 max(time)+100])
    xlim([max(time)-yr_1800_CE max(time)])                   

    h19a=gca; 
    set(h19a,'FontSize',10)

    h19c = gcf;
    fig19name = [params.out_name, '_smoothed_NPP_NCB_H_M_final'];
    saveas(h19c, fig19name,'jpg');

end


%------------------------------------
% Diagnostics for soil temperature
% -----------------------------------

% ALT:
figure(33)
subplot(4,1,1)
plot(time,-ann_ALD1_max, time,-ann_ALD2_max, time,-ann_ALD3_max, time,-ann_ALD_max)
xlabel('\fontsize{10}Time')
title('\fontsize{10}Active Layer Depth comparisons');
legend('\fontsize{10}ALD1 @ 0','\fontsize{10}ALD2 @ -0.08', '\fontsize{10}ALD3@ -0.16',...
    '\fontsize{10}Final ALD')
legend('orientation','Horizontal','Location','South')


layers_of_int = [1, 5, 11, 21, 41];
subplot(4,1,2)
plot(soil_node_temp_month_saveMAT(1:(100*12),layers_of_int))
title('\fontsize{10}Soil Temps, first 100 years of simulation');
legend('\fontsize{10}0 cm','\fontsize{10}20 cm', '\fontsize{10}50 cm',...
    '\fontsize{10}100 cm', '\fontsize{10}300 cm')
legend('orientation','Horizontal','Location','South')
subplot(4,1,3)
plot(soil_node_temp_month_saveMAT(((sim_start+55)*12):((sim_start+65)*12),layers_of_int)); hold on
plot(mon_snowdepth(((sim_start+55)*12):((sim_start+65)*12))*100); 
plot(-mon_wtd(((sim_start+55)*12):((sim_start+65)*12))*100); hold off
title('\fontsize{10}Monthly soil temps, 2005-2015');
legend('\fontsize{10}0 cm','\fontsize{10}20 cm', '\fontsize{10}50 cm',...
    '\fontsize{10}100 cm', '\fontsize{10}300 cm', '\fontsize{10}SnowDepth_cm', '\fontsize{10}WTD')
legend('Location','East')
subplot(4,1,4)
plot(ann_Z_total); hold on
plot(ann_Z_total - ann_ALD_max)
plot(ann_Z_total - mon_wtd(jul))
title('\fontsize{10}Peat Height, WT, ALD');
legend('\fontsize{10}Peat height','\fontsize{10}Annual ALD', '\fontsize{10}July WTD')
legend('orientation','Horizontal','Location','South')
hold off

h33c = gcf;
fig33name = [params.out_name, '_soilT_metrics'];
saveas(h33c, fig33name,'jpg');

