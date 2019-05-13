function [mon_runon] = hpm20_mon_runon(WTD, THICK, DEPTH, layer_frozen, params)

% monthly runon in m water depth

% % ******************
% *** RUN-ON ***
% ******************

% function calculates monthly water runon (m/y) (just a placeholder for now)
% RUNON PARAMETERS (values read in params script)
%  ? runon_c1: total peat height where runon declines by 50% (set deep for fen)
%  ? runon_c2: parameter for rate of decline of runon (see ?HPM vegetation productivity.xls?)
%  ? runon_c3: magnitude of maximum annual runon (m)
%      modify for monthly?

peatdepth = sum(THICK);

if (params.runon_c3 > 0)
    Ron = (params.runon_c3) * (1 - 0.5 * ( 1 + erf((peatdepth - params.runon_c1)/(sqrt(2)*params.runon_c2))));
end

% TESTING A ?NEW(?)? IDEA (APRIL 2009)
% run-on declines with decreasing WTD (i.e., as water table rises) 
% Ron factor rises from zero at WTD = 0 to one at WTD = 0.1 m
% should this modify or replace runon calculation of previous few lines?

if (params.runon_c3 > 0)
    runon_c4 = 0.0; % WTD depth above which run-on is zero
    runon_c5 = 0.1; % WTD depth below which run-on is maximum
    Ron = Ron * min(1.0, max(0.0, (WTD-runon_c4)/runon_c5));  
else
    Ron = 0;
end

mon_runon = Ron;

% END RUNON FUNCTION

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


