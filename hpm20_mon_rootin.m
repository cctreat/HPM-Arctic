function root_in = hpm20_mon_rootin(depthvec, thickvec, params, nppvec, zwt, alt, peatheight, onevct)
% root_in = hpm20_mon_rootin(depthvec, thickvec, params, nppvec, zwt, alt, peatheight, onevct)

% to determine root inputs by PFT for all cohorts in the root zone
% root profiles follow ideas of Bauer (2004) 
% zwt and alt are growing season maximum monthly values

sedge_tot_root = params.bg_frac_npp .* nppvec .* params.sedges; %
non_sedge_tot_root = params.bg_frac_npp .* nppvec .* (params.vasculars - params.sedges);  

% FOR SENSITIVITY TEST: make all roots 'sedge'
% sedge_tot_root = params.bg_frac_npp .* params.vasculars;
% non_sedge_tot_root = params.bg_frac_npp .* 0;

% FOR SENSITIVITY TEST: make all roots 'non-sedge'
% non_sedge_tot_root = params.bg_frac_npp .* params.vasculars;
% sedge_tot_root = params.bg_frac_npp .* 0;

% determine root input depth ?zstar?
z1 = max(zwt, params.rootin_min);  % root input to WT or minimum depth
z2 = params.rootin_max;   % max root input depth if deep WT; then restricted to active layer

if (params.pf_flag > 0.5) % permafrost
    z2 = min(alt, params.rootin_max);
end
zstar = min(z1, z2);

if (params.pf_flag > 0.5)  % permafrost
    upper_sedge_rootin_frac = min(1., peatheight/ (alt + eps));  % fraction sedge shallow root input into peat, if shallower than active layer
    non_sedge_rootin_frac = min(1., peatheight/ (alt + eps));  % fraction root input into peat, if shallower than active layer
else   % no permafrost
    upper_sedge_rootin_frac = min(1., peatheight/ params.rootin_d80);  % fraction sedge shallow root input into peat, if shallower than rooting
    non_sedge_rootin_frac = min(1., peatheight / (zstar + eps));  % fraction root input into peat, if shallower than rooting
end

% SF: new routines for root input (August 2011) uniform input into cohorts, rather than by thickness ? WHY???
%  from Julie Talbot?s notes: 
%      There was runaway peat accumulation when we did simulations without mosses for the Turetsky paper. 
%      We thought it might be related to the root input. 
%      The changes you made to the code then (in August 2011) fixed the problem. 

% ***SEDGE ROOTS***
%  sedge roots: uniform input per layer for upper range of sedge roots (depth < 'd80' parameter (depth to 80% of root input)
%  sedge roots: input proportional to layer thickness below 'd80', with total of 20% from 'd80' to 2 meters (params.rootin_c4)

input_equal_per_layer = 1;

if (input_equal_per_layer > 0.5)
    
    number_root_layers = find(depthvec > params.rootin_d80, 1,'first')-1;
 
    if (isempty(number_root_layers))
         number_root_layers = find(thickvec > 0, 1,'last');
         if (isempty(number_root_layers))
             number_root_layers = 1;
         end
    end

    upper_sedge_root_in = 0.8 / (max(1,number_root_layers));
    tf_root1 = thickvec > 0;
    tf_root2 = depthvec <= params.rootin_d80;
    tf_root3 = tf_root1 .* tf_root2;
    if (isempty(tf_root3))
        tf_root3 = 0 * thickvec;
        tf_root3(1) = 1;
    end
    upper_sedge_root_in = upper_sedge_root_in .* tf_root3 * upper_sedge_rootin_frac;

    lower_sedge_root_in = (1 - 0.8)/(exp(-params.rootin_d80) - exp(-params.rootin_sedge_max)) * exp(-depthvec) .* thickvec;
%     * equation above has exp. decay of input with total sum to bottom (rootin_c4) of 0.2 (numerator)
    tf_root4 = depthvec <= 2;  % no sedge roots below 2 meters
    lower_sedge_root_in = lower_sedge_root_in .* tf_root4;
    
    sedge_root_in = upper_sedge_root_in + lower_sedge_root_in;
    
else   %%% original algorithm
    
    sedge_root_in = thickvec .* (params.rootin_alpha * exp(-params.rootin_alpha * depthvec))...
                .* (onevct - min(onevct,fix(depthvec/params.rootin_c4)));

    norm1 = sum(sedge_root_in) / min(1,(1-exp(-params.rootin_alpha * peatheight)));
            
% sedge_root_in = sedge_root_in / (sum(sedge_root_in) + eps);  % normalize total to 1.0??
    sedge_root_in = sedge_root_in / (norm1 + eps);  % normalize total to 1.0
  
end

% -- if permafrost, compress sedge root input to active layer

if (params.pf_flag > 0.5)
    zbottom = cumsum(thickvec);
    tf_zbot = zbottom < alt;
    sedge_root_in_pf = sedge_root_in .* tf_zbot; % restrict root input to Active Layer
    if (sum(sedge_root_in_pf) > 0)
        pf_factor = 1/sum(sedge_root_in_pf);
    else
        pf_factor = 1.;
    end
    sedge_root_in_pf = sedge_root_in_pf * pf_factor;
    sedge_root_in = sedge_root_in_pf;
end

% ***NON-SEDGE ROOTS***
%  non-sedge roots: uniform input per layer for (rather than proportional to layer thickness)

if (input_equal_per_layer > 0.5)
    
    number_root_layers = find(depthvec > zstar, 1,'first') - 1;

    if (isempty(number_root_layers))
         number_root_layers = find(thickvec > 0, 1,'last');
         if (isempty(number_root_layers))
             number_root_layers = 1;
         end
    end

    non_sedge_root_in = onevct / (max(1,number_root_layers));
%    non_sedge_root_in = root_frac * (onevct / (number_root_layers));
    tf_root1 = thickvec > 0;
    tf_root2 = depthvec <= zstar;
    tf_root = tf_root1 .* tf_root2;
    if (isempty(tf_root))
        tf_root = 0 * thickvec;
        tf_root(1) = 1;
    end
    non_sedge_root_in = non_sedge_root_in .* tf_root * non_sedge_rootin_frac;

else
    
% original version (below) uses error function to get a smooth boundary,second has uniform input to zstar
% second version lost about 5% of root mass due to discretization(?), hence divided by sum...

%     non_sedge_root_in = (thickvec/zstar) .* (onevct - 0.5*(onevct + erf((depthvec - zstar*onevct)/(sqrt(2)*params.rootin_c5))));

    non_sedge_root_in = (thickvec/zstar) .* (depthvec < zstar);

    non_sedge_root_in = non_sedge_root_in / (sum((thickvec/zstar) .* (depthvec < zstar)) + eps);  % normalize total to 1.0??
    non_sedge_root_in = non_sedge_root_in * min(1, peatheight/zstar);  % adjust total to fraction of root zone that is peat
    
%     norm2 = sum((thickvec/zstar) .* (depthvec < zstar)) / min(1,(peatheight/zstar));
%     non_sedge_root_in = non_sedge_root_in / (norm2 + eps);  % normalize total to 1.0

end

root_in = sedge_root_in * sedge_tot_root + non_sedge_root_in * non_sedge_tot_root;

%j5a = [ norm1 norm2 ]

return
