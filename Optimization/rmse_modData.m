% Function to minimize
function rmseOut = rmse_modData(modIn, dataIn)
% function for root mean square error

rmseOut = sqrt(sum((modIn - dataIn).^2));
return