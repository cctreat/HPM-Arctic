function [mon_runoff transmis] = hpm20_mon_runoff(WTD, THICK, DEPTH, DENS, ONEVEC, layer_frozen, params, ALT)
% function calculates monthly runoff in m water depth
 
% ------------

% first calculate hydraulic conductivity of cohorts down profile
%    hydraulic conductivity as function of bulk density from Radforth (1977)

% Calculate hydraulic conductivity down the cohort profile
% log_10(K) = (150 - 3*rho)/70    % From Radforth (1977)
% 150/70 = 2.142857 ; 3/70 = 0.042857; 10^(150/70) = 138.94955 ; 10^(-3/70) = 0.906031

Log10HydrCond = 2.14287 * ONEVEC - 0.042857 * DENS;
HydrCond = exp(log(10)*Log10HydrCond);

% ALTERNATIVE CALCULATION (without log and exp)
% HydrCond = 138.94955 * (ONEVEC * 0.906031) .^ DENSITY;

% ------------

% second calculate location of WT and AL 
%        transmissivity a function of thickness of saturated zone (above pf if present)

trans_counter_wt = find(DEPTH > WTD,1);
if (isempty(trans_counter_wt))
    trans_counter_wt = find(DEPTH > 0,1,'last');
end

if (params.pf_flag < 0.5)   % no permafrost

    transmis = sum(THICK(trans_counter_wt:end) .* HydrCond(trans_counter_wt:end)) / ...
                   sum(THICK(1:end) .* HydrCond(1:end));    

else   %  permafrost

    trans_counter_pf = find(DEPTH > ALT,1);
    if (isempty(trans_counter_pf))
        trans_counter_pf = find(DEPTH < ALT,1,'last');
        trans_counter_pf = max(trans_counter_pf,1);
    end

    trans_counter_wt1 = min(trans_counter_wt, trans_counter_pf);

    transmis = sum(THICK(trans_counter_wt1:trans_counter_pf) .* HydrCond(trans_counter_wt1:trans_counter_pf)) / ...
                   sum(THICK(1:trans_counter_pf) .* HydrCond(1:trans_counter_pf));   

end

transmis = transmis * (sum(layer_frozen(1:3)) < 3);  % no runoff is surface peat is frozen?

% Annual version scaled transmissivity range from 0.5 to 1.0
%  ***dropping this from monthly time step

% transmis = params.Roff_c3 + (1 - params.Roff_c3) * transmis;  

% ALTERNATIVE METHOD #2
% below = DEPTH > WTD;  % creates 'logical' vector with ones for elements with depth>WTD, zeros otherwise
% transmissivity = sum(THICK .* (HydrCond .* below)) / sum(THICK .* HydrCond);

% ALTERNATIVE METHOD #1
% NOTE: this first one gives different results--why?
% THICK1 = THICK;
% THICK1(THICK1<=WTD) = 0;
% transmissivity = sum(THICK1 .* HydrCond) / sum(THICK .* HydrCond);

% ------------

% function calculates annual water runoff (m/y)
% follows revised version of PAM from Roulet (see ref above)
% reference for PAM is:
 % Hilbert D, NT Roulet, TR Moore. 2000. Modelling and analysis of peatlands as dynamical systems, J. Ecol. 88:230-242.


% RUNOFF PARAMETERS (values read in params script)
%  ? Roff_c1: target annual runoff (= annual precip minus annual ET plus about 10 cm/y)
%      modify for monthly?
%  ?  Roff_c2: linear increase in runoff (m/y) per meter of total peat height)
%  ?  Roff_c2a: offset in total peat height to start peat accumulation in a low area.
%  ?  Roff_c3: minimum profile relative transmissivity (set to 0.5 for annual simulation).
%      modify for monthly?
%  ?  Roff_c4: threshold WTD (negative, so inundation) for immediate spillover (= Roff_c4 ? WTD); now used in main code.

months_with_runoff = 12;

runoff1 = params.Roff_c1/months_with_runoff * (1 + params.Roff_c2 * (sum(THICK) - params.Roff_c2a));  % modified run-off (March 2010)
Roff = transmis * runoff1;

%Roff = Roff * 0.75;
% eliminate the possibility of negative runoff
if Roff < 0
    Roff = 0;
end

mon_runoff = Roff;

% END RUNOFF FUNCTION

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


