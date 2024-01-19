%Given allocatinon only in env A
%Mirror image fungi, but varying differences between them 


% set parameter values randomly 
% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

rtot = .8; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 8000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 
results3 = results2B; 

env_period = 365; 

figure 
clf 

    difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

preference=1; 

d1_1 = .1;
d2_2 = .1;

count = 1;
for init = 1:2

    if init == 1
        x0(3) = 1; %F1
    else
        x0(3) = 10; %F1
        tspan = [1 80000]; 
    end

    for j = 1:3

        if j == 1
            d2_1 = .05;
            d1_2 = .05;
        elseif j == 2
            d2_1 = .1;
            d1_2 = .1;
        else
            d2_1 = .15;
            d1_2 = .15;
        end


        for i = 1:length(propA_vals)
            propA = propA_vals(i);
            envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

            %first run simulation for reward strategy 100% fungus 1 in both environments
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*preference, rtot*preference, rtot*(1-preference), rtot*(1-preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(1,i) = mean(final_res(1,:));
            results2A(1,i) = median(final_res(1,:));
            results2A(2,i) = quantile(final_res(1,:), .25);
            results2A(3,i) = quantile(final_res(1,:), .75);


            %next run simulation for reward strategy 100% fungus 2 in both
            %environments
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(1-preference), rtot*(1-preference), rtot*preference, rtot*preference, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(2,i) = mean(final_res(1,:));

            results2B(1,i) = median(final_res(1,:));
            results2B(2,i) = quantile(final_res(1,:), .25);
            results2B(3,i) = quantile(final_res(1,:), .75);

            %run simulation for reward strategy 50:50
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            results(3,i) = mean(final_res(1,:));

            results3(1,i) = median(final_res(1,:));
            results3(2,i) = quantile(final_res(1,:), .25);
            results3(3,i) = quantile(final_res(1,:), .75);

        end

        subplot(2,3,count)

        hold on
        boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
        boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
        boundedline(propA_vals, results3(1,:), [results3(1,:)-results3(2,:);results3(3,:)- results3(1,:)]', 'linewidth', 2, 'color', 'k', 'alpha', 'transparency', 0.5)

        xlabel('Proportion of time in environment A')
        
        count = count+1;

    end
end


subplot(2,3,1)
title('d_{i,j} < d_{i,i}')

subplot(2,3,2)
title('d_{i,j} = d_{i,i}')
legend({'Bet hedging'; ''; 'F_1 preference'; ''; 'F_2 preference'; ''; ''; ''})

subplot(2,3,3)
title('d_{i,j} > d_{i,i}')
