% Optimization



%1. Run model for soil temperature routine
hpm20_mon_params_Selwyn_win.m
hpm20_mon_main.m

%2. clear
clear all

%3. Load optimization parameters
hpm20_mon_params_Selwyn_winOptim

%4. run optimization

%5. Update parameterization & re-run models to update soil temperature
%profile
hpm20_mon_params_Selwyn_win.m
hpm20_mon_main.m


%6. clear
clear all

%7. Load optimization parameters
hpm20_mon_params_Selwyn_winOptim

%8. re-run optimization