%Given allocatinon only in env A
%Mirror image fungi, but varying differences between them 


% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .05; %allocation of Carbon to Carbon pool
s = 0.001; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 3000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 

env_period = 365; 

for d = 1:3; 
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 1, 1, 0, 0, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 

    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 0, 0, 1, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    if propA == .5; 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
            title('100% Fungus 2 preference')

            xticks([0:365:3000])
            xlabel('Days')
    end

     % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

subplot(2,3,d)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 18])
legend({'100%A'; '100%B'; '50:50'})
end

%% Try again but with median and std rather than simulation over time

%Mirror image fungi, but varying differences between them 


% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .05; %allocation of Carbon to Carbon pool
s = 0.001; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 3000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2_50 = nan(3, length(propA_vals)); 


env_period = 365; 

for d = 1:3; 
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 1, 1, 0, 0, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 
    results2A(1,i) = median(final_res(1,:)); 
    results2A(2,i) = quantile(final_res(1,:), .25); 
    results2A(3,i) = quantile(final_res(1,:), .75); 


    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 0, 0, 1, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

     % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

    results2_50(1,i) = median(final_res(1,:)); 
    results2_50(2,i) = quantile(final_res(1,:), .25); 
    results2_50(3,i) = quantile(final_res(1,:), .75); 

end

subplot(2,3,d)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 18])
legend({'100%A'; '100%B'; '50:50'})

subplot(2,3,d+3)
plot(propA_vals, results2A(1,:))
hold on 
plot(propA_vals, results2A(2,:), 'color', [0 0.4470 0.7410], 'linestyle', '--')
plot(propA_vals, results2A(3,:), 'color', [0 0.4470 0.7410], 'linestyle', '--')
plot(propA_vals, results2_50(1,:), 'color', [0.9290 0.6940 0.1250])
hold on 
plot(propA_vals, results2_50(2,:), 'color', [0.9290 0.6940 0.1250], 'linestyle', '--')
plot(propA_vals, results2_50(3,:), 'color', [0.9290 0.6940 0.1250], 'linestyle', '--')

ylabel('Tree biomass')
xlabel('Proportion of time in environment A')
ylim([0 18])


legend({'100% A'; ''; ''; '50:50'})


end


%% Leakiness in different environments 
% now vary reward strategy in fixed environmental regime

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .05; %allocation of Carbon to Carbon pool
s = 0.001; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
leakiness_vals = 0:.05:1; 
env_period_vals = [180 365 365*2 365*3];
propA = .5; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

tspan = [1 8000]; 

for d = 1:3; 
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1-leakiness; 

    for j = 1:length(env_period_vals)

        env_period = env_period_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));
    end


end

figure(1)
subplot(1,3,d)
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend({'6 months'; '1 year'; '2 years'; '3 years'}); 
title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 16])

end

%%


