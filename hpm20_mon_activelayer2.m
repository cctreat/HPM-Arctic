function [ALD1, ALD2, ALD3] = hpm20_mon_activelayer(soilT_yearMoxNode, params_gipl)

% COMPUTE ALDs FROM MEAN MONTHLY SOIL NODE TEMPS
% Using technical definition of permafrost: layer frozen for last 2 years.
% ALD1 is warm end of freezing zone band
% ALD2 is 'middle' of freezing zone band
% ALD3 is cold end of freezing zone band

%soilT_yearMoxNode = reshape(soilT_yearMoxNode, [36, 63]); % reshape soil temp array
maxNodeT_2yr = max(soilT_yearMoxNode, [], 1);                % find maximum soil temperature over the last 2 years

TF_nodeFroz_2yr_1 = maxNodeT_2yr < (params_gipl.Tfr + params_gipl.FIT); % compare with freezing metric 1 (0)
if sum(TF_nodeFroz_2yr_1) > 0
    ALD1 = params_gipl.soilNodeDepth(find(TF_nodeFroz_2yr_1, 1));
else ALD1 = 10;
end
   
TF_nodeFroz_2yr_2 = maxNodeT_2yr < (params_gipl.Tfr); % middle of thawing zone (above bedrock) 
if sum(TF_nodeFroz_2yr_2) > 0
    ALD2 = params_gipl.soilNodeDepth(find(TF_nodeFroz_2yr_2, 1));
else ALD2 = 10;
end

TF_nodeFroz_2yr_3 = maxNodeT_2yr < (params_gipl.Tfr - params_gipl.FIT); % Bottom of thawing zone
if sum(TF_nodeFroz_2yr_3) > 0
    ALD3 = params_gipl.soilNodeDepth(find(TF_nodeFroz_2yr_3, 1)) ;
else ALD3 = 10;
end


return

