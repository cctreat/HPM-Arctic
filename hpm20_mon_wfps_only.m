function [new_wfps, new_PEAT_wat] = ...
                  hpm20_mon_wfps_only(WTD, Zstar,THICK,DEPTH,POROSITY, ONEVEC,ZEROVEC,params ,ALT)

% Computes peat WFPS only, if WTD is already known or prescribed

% OUTPUT VARIABLES
% ------------------
% new_wfps is updated peat cohort WFPS profile 
% new_PEAT_wat is final water content of peat profile (saturated + unsaturated) [m or m3/m2]

% INPUT VARIABLES
% ------------------
% WTD is given water table depth below surface [m]
% Zstar is factor for unsaturated WFPS, and is function of bulk density
% DENS is cohort bulk density profile [kg/m3]
% THICK is vector of cohort thicknesses (m)
% DEPTH is vector of cohort mid-point depths (m, positive down)
% POROSITY is vector of cohort porosities (- or m3/m3)
% ONEVEC & ZEROVEC are vectors of ones and zeros
% params is model input parameters
% ALT is current active layer thickness (is this needed?)
% ------------------

% COMPUTED VARIABLES
% ------------------
% zwtd is cohort distance above water table (used to compute WFPS of unsaturated peat) 
% ------------

% compute peat water content for this given WT position estimate

if (WTD <= 0.)
    new_wfps = ONEVEC;
    new_PEAT_wat = sum(POROSITY .* THICK);

else

    zwtd = DEPTH - WTD;
    zwtd = max(ZEROVEC, -zwtd);
    new_wfps = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_PEAT_wat = sum(new_wfps .* POROSITY .* THICK);

end

return;   %  ADD any checks on WT position (e.g., is it within the peat?)
% ------------
