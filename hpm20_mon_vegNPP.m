function npptotalmax = hpm20_mon_vegNPP(wtd_opt,wtd_range,pd_opt,pd_range,npp_rel,k0,num_veg)
% function tot_npp = hpm20_mon_vegNPP (wtd_opt,wtd_range,pd_opt,pd_range,npp_rel,k0,num_veg)

% function is called from the parameters script, so cannot use the ?params? structure

% v.20: making version for version control - June 2015

% v12 for arbitary number of PFTs

% v9

% v8 no change from v6

% function generates plots of veg NPP and calculates total NPP as function of 2-dimensional ?control? space

% each vegetation type will have an optimal position, and a surface of diminishing NPP

% response modeled as 2-D Gaussian functions, with (in most cases)
%    different variances on either side of the optimum to skew the function
%    to have a long (productivity) tail in one direction (e.g., increasing WTD)
%    and a short (productivity) tail in the other direction (e.g., decreasing WTD)

% PARAMETERS
%-----------

%  WTD_opt = optimum WTD
%  a1 = NPP sensitivity to WTD increasing (deeper WT) (Gaussian curve variance)
%  a2 = NPP sensitivity to WTD decreasing (shallower WT) (Gaussian curve variance)

%  PD_opt = optimum PD or ALD
%  b1 = NPP sensitivity to PD increasing (deeper peat profile) (Gaussian curve variance)
%  b2 = NPP sensitivity to PD decreasing (Gaussian curve variance)

%  npp_rel = relative maximum NPP from input

[X,Y] = meshgrid(-.5:0.1:1, 0:0.1:8);
minarray = 0.00001*ones(size(X));
NPPmax = zeros(13,1);

XX = [-0.5:0.1:1];
YY = [0:0.1:8];

% loop through PFTs

for nveg = 1:1:num_veg

    atest1 = XX < wtd_opt(nveg);  % one if WTD <  optimum, zero otherwise
    atest2 = XX >= wtd_opt(nveg); % one if WTD >= optimum, zero otherwise
    btest1 = YY < pd_opt(nveg);   % one if PD  <  optimum, zero otherwise
    btest2 = YY >= pd_opt(nveg);  % one if PD  >= optimum, zero otherwise
    btest1 = btest1';
    btest2 = btest2';
    
    NPPa1 = exp(-(((XX - wtd_opt(nveg)) / wtd_range(1,nveg)).^2));
    NPPa2 = exp(-(((XX - wtd_opt(nveg)) / wtd_range(2,nveg)).^2));
    NPPb1 = exp(-((YY - pd_opt(nveg)) / pd_range(1,nveg)).^2);
    NPPb2 = exp(-((YY - pd_opt(nveg)) / pd_range(2,nveg)).^2);
    NPPb1 = NPPb1';
    NPPb2 = NPPb2';
    
    NPP1 = npp_rel(nveg) * ((NPPb1 .* btest1 + NPPb2 .* btest2) * (NPPa1 .* atest1 + NPPa2 .* atest2));
    NPPmax(nveg) = max(max(NPP1));
    NPP(nveg,:,:) = max(NPP1, minarray);
    
    if (nveg == 1)
        NPPtotal = NPP(1,:,:);
    else
        NPPtotal = NPPtotal + NPP(nveg,:,:);
    end
end

NPPtotal = squeeze(NPPtotal);
NPPtotal = max(NPPtotal, minarray);
npptotalmax = max(max(NPPtotal));

% TOTAL of mosses, sedges, vasculars, trees, and peat (decomp weighted)
%  can be computed for evaluation, but not returned to main code

% for nveg = 1:1:num_veg
% 
%     if (nveg == 1) 
%         NPPmosstotal = NPP(nveg,:,:) * params.mosses(nveg);
%         NPPsedgetotal = NPP(nveg,:,:) * params.sedges(nveg);
%         NPPtreetotal = NPP(nveg,:,:) * params.trees(nveg);
%         NPPvasctotal = NPP(nveg,:,:) * params.vasculars(nveg);
%         NPPpeattotal = NPP(nveg,:,:) / (k0(nveg) + eps);
%     else
%         NPPmosstotal = NPPmosstotal + NPP(nveg,:,:) * params.mosses(nveg);
%         NPPsedgetotal = NPPsedgetotal + NPP(nveg,:,:) * params.sedges(nveg);
%         NPPtreetotal = NPPtreetotal + NPP(nveg,:,:) * params.trees(nveg);
%         NPPvasctotal = NPPvasctotal + NPP(nveg,:,:) * params.vasculars(nveg);
%         NPPpeattotal = NPPpeattotal + NPP(nveg,:,:) / (k0(nveg) + eps);
%     end
% 
% end
% 
% NPPpeattotal = squeeze(NPPpeattotal);
% NPPmosstotal = squeeze(NPPmosstotal);
% NPPvasctotal = squeeze(NPPvasctotal);
% NPPtreestotal = squeeze (NPPtreestotal);
% 
% NPPpeat = max(NPPpeat, minarray);
% NPPmosstotal = max(NPPmosstotal, minarray);
% NPPvasctotal = max(NPPvasctotal, minarray);
% NPPtreestotal = max (NPPtreestotal, minarray);
% 
% % 
% npppeattotalmax = max(max(NPPpeattotal));
% nppmossmax = max(max(NPPmosstotal))
% nppvascmax = max(max(NPPvasctotal))
% npptreesmax = max(max(NPPtreestotal))
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot results

% plot or not

plot_flag = 1;

if (plot_flag > 0)

  % ---------  
  % figure(21)
  % ---------  

    figure(21)
    
    for nveg = 1:1:num_veg
        subplot(2,3,nveg)
        contourf(X, Y, squeeze(NPP(nveg,:,:)),100,'LineStyle','none')
        xlim([-0.1 0.6])

%         str = sprintf('PFT #', nveg);
%         title(str, 'FontSize',14)   
        title(['\fontsize{14}PFT #',num2str(nveg)]);

        ylim([0 5])
       zlim([0 1])
        colorbar
        caxis([0 2])
        xlabel('\fontsize{14}water table depth [m]')
        ylabel('\fontsize{14}peat or AL depth [m]')
        zlabel('\fontsize{14}relative NPP')
        h1a=gca;
        set(h1a,'FontSize',14)
    end
    
    subplot(2,3,nveg+1)
    contourf(X, Y, NPPtotal,100,'LineStyle','none')
    xlim([-0.1 0.6])
    title('\fontsize{14}all PFT total');
    ylim([0 5])
    zlim([0 1])
    colorbar
    caxis([0 2])
    xlabel('\fontsize{14}water table depth [m]')
    ylabel('\fontsize{14}peat or AL depth [m]')
    zlabel('\fontsize{14}relative NPP')
    h1a=gca; 
    set(h1a,'FontSize',14)

end

return

