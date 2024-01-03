%%

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.01;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.01;  %density dependent mortality effect of F2 on F2
rtot = .8;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

env_period_vals = [4*365];


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
propA_vals = 0:.02:1; 
results = nan(3,length(propA_vals)); 

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, rtot, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
       final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 

    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, 0, rtot, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
      final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l,  rtot/5, rtot/5,  rtot/5,  rtot/5, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
      final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

figure
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')


%%
%let's look across environemtal regimes and see how "optimal" leakiness
%affects fungal and host biomass dynamics compared to .. 0% leakiness? 

% Response Strategy (Whichever Fungus is best right now) 
% Loyal Strategy (Fungus A Only) 
% Optimally Leaky (A little of both to maximize average biomass) 

% Want biomass timeseries for these three strategies across PropA values, 

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.04;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.04;  %density dependent mortality effect of F2 on F2
rtot = .2;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

env_period_vals = [4*365];

results = nan(length(env_period_vals), length(leakiness_vals)); 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

propA_vals = .8 %[0 .2 .4 .6 .8 1]; 
for p = 1:length(propA_vals)
   propA = propA_vals(p);

   env_period = env_period_vals(j);
   envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

%Find optimal leakiness
optimal_leakinesss_sol = []; 
optimal_leakiness_biomass = 0; 
optimal_leakiness_value = []; 

leakiness_vals = 0:.05:1; 
results = nan(1, length(leakiness_vals));
for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i) = mean(final_res(1,:)); 

        if optimal_leakiness_biomass < mean(final_res(1,:))
            optimal_leakiness_biomass = mean(final_res(1,:)); 
            optimal_leakiness_sol = final_res; 
            optimal_leakiness_value = leakiness; 
        end

end


% Now Responsive solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   responsive_sol = deval(sol, tspan(2)-env_period*3:tspan(2));
   

% Now Loyal solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, rtot, 0, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   loyal_sol = deval(sol, tspan(2)-env_period*3:tspan(2));

end

figure(2)
plot(leakiness_vals, results)
xlabel('Leakiness')
ylabel('Average tree biomass')

figure(1)
subplot(1,3,1)
b = plot(tspan(2)-env_period*3:tspan(2), loyal_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title('Loyal')

subplot(1,3,2)
b = plot(tspan(2)-env_period*3:tspan(2), responsive_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title('Responsive')

subplot(1,3,3)
b = plot(tspan(2)-env_period*3:tspan(2), optimal_leakiness_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title(['Optimally Leaky (' num2str(optimal_leakiness_value.*100) '%)'])

%%



% Responsive Strategy (Whichever Fungus is best right now) 
% Loyal Strategy (Fungus A Only) 
% Optimally Leaky (A little of both to maximize average biomass) 

% Want biomass timeseries for these three strategies across PropA values, 

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.02;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.02;  %density dependent mortality effect of F2 on F2
rtot = .2;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

env_period_vals = [4*365];

results = nan(length(env_period_vals), length(leakiness_vals)); 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

propA_vals = 0:.05:1; 
propA_results = nan(12, length(propA_vals)); 
optimal_leakiness_values = nan(1, length(propA_vals)); 

for p = 1:length(propA_vals)
   propA = propA_vals(p);

   env_period = env_period_vals(j);
   envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

%Find optimal leakiness
optimal_leakinesss_sol = []; 
optimal_leakiness_biomass = 0; 
optimal_leakiness_value = []; 

leakiness_vals = 0:.05:1; 
results = nan(1, length(leakiness_vals));
for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i) = mean(final_res(1,:)); 

        if optimal_leakiness_biomass < mean(final_res(1,:))
            optimal_leakiness_biomass = mean(final_res(1,:)); 
            optimal_leakiness_sol = final_res; 
            optimal_leakiness_value = leakiness; 
        end

end

optimal_leakiness_values(p) = optimal_leakiness_value; 

% Now Responsive solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   responsive_sol = deval(sol, tspan(2)-env_period*3:tspan(2));
   

% Now Loyal solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, rtot, 0, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   loyal_sol = deval(sol, tspan(2)-env_period*3:tspan(2));



propA_results(1, p) = mean(optimal_leakiness_sol(1,:)); 
propA_results(2,p) = std(optimal_leakiness_sol(1,:)); 
propA_results(3,p) = mean(optimal_leakiness_sol(3,:)); 
propA_results(4,p) = std(optimal_leakiness_sol(3,:)); 

propA_results(5, p) = mean(responsive_sol(1,:)); 
propA_results(6,p) = std(responsive_sol(1,:)); 
propA_results(7,p) = mean(responsive_sol(3,:)); 
propA_results(8,p) = std(responsive_sol(3,:)); 

propA_results(9, p) = mean(loyal_sol(1,:)); 
propA_results(10,p) = std(loyal_sol(1,:)); 
propA_results(11,p) = mean(loyal_sol(3,:)); 
propA_results(12,p) = std(loyal_sol(3,:)); 

end

figure(3) 
subplot(1,2,1)
plot(propA_vals, propA_results(1,:), 'color', [0 1 0]) 
hold on 
plot(propA_vals, propA_results(1,:)+ propA_results(2,:), 'linestyle', '--', 'color', [0 1 0]) 
plot(propA_vals, propA_results(1,:)- propA_results(2,:), 'linestyle', '--', 'color', [0 1 0]) 

plot(propA_vals, propA_results(5,:), 'color', [1 0 0]) 
hold on 
plot(propA_vals, propA_results(5,:)+ propA_results(6,:), 'linestyle', '--', 'color', [1 0 0]) 
plot(propA_vals, propA_results(5,:)- propA_results(6,:), 'linestyle', '--', 'color', [1 0 0]) 

plot(propA_vals, propA_results(9,:), 'color', [0 0 1]) 
hold on 
plot(propA_vals, propA_results(9,:)+ propA_results(10,:), 'linestyle', '--', 'color', [0 0 1]) 
plot(propA_vals, propA_results(9,:)- propA_results(10,:), 'linestyle', '--', 'color', [0 0 1]) 

ylabel('Tree Biomass')
xlabel('Proportion of time in environment A')

subplot(1,2,2)
plot(propA_vals, propA_results(3,:), 'color', [0 1 0]) 
hold on 
plot(propA_vals, propA_results(3,:)+ propA_results(4,:), 'linestyle', '--', 'color', [0 1 0]) 
plot(propA_vals, propA_results(3,:)- propA_results(4,:), 'linestyle', '--', 'color', [0 1 0]) 

plot(propA_vals, propA_results(7,:), 'color', [1 0 0]) 
hold on 
plot(propA_vals, propA_results(7,:)+ propA_results(8,:), 'linestyle', '--', 'color', [1 0 0]) 
plot(propA_vals, propA_results(7,:)- propA_results(8,:), 'linestyle', '--', 'color', [1 0 0]) 

plot(propA_vals, propA_results(11,:), 'color', [0 0 1]) 
hold on 
plot(propA_vals, propA_results(11,:)+ propA_results(12,:), 'linestyle', '--', 'color', [0 0 1]) 
plot(propA_vals, propA_results(11,:)- propA_results(12,:), 'linestyle', '--', 'color', [0 0 1]) 

ylabel('Fungus 1 Biomass')
xlabel('Proportion of time in environment A')

figure(4)
plot(propA_vals, optimal_leakiness_values); 
ylabel('Optimal Leakiness')
xlabel('Proportion of time in environment A')

% figure(2)
% plot(leakiness_vals, results)
% xlabel('Leakiness')
% ylabel('Average tree biomass')
% 
% figure(1)
% subplot(1,3,1)
% b = plot(tspan(2)-env_period*3:tspan(2), loyal_sol);
% b(1).Color = [0 .9 .3];
% b(2).Color = [0 .8 .7];
% b(3).Color = 'r';
% b(4).Color = 'b';
% b(5).Color = [.5 0 .5];
% 
% b(2).LineStyle = '--';
% b(5).LineStyle = '--';
% title('Loyal')
% 
% subplot(1,3,2)
% b = plot(tspan(2)-env_period*3:tspan(2), responsive_sol);
% b(1).Color = [0 .9 .3];
% b(2).Color = [0 .8 .7];
% b(3).Color = 'r';
% b(4).Color = 'b';
% b(5).Color = [.5 0 .5];
% 
% b(2).LineStyle = '--';
% b(5).LineStyle = '--';
% title('Responsive')
% 
% subplot(1,3,3)
% b = plot(tspan(2)-env_period*3:tspan(2), optimal_leakiness_sol);
% b(1).Color = [0 .9 .3];
% b(2).Color = [0 .8 .7];
% b(3).Color = 'r';
% b(4).Color = 'b';
% b(5).Color = [.5 0 .5];
% 
% b(2).LineStyle = '--';
% b(5).LineStyle = '--';
% title(['Optimally Leaky (' num2str(optimal_leakiness_value.*100) '%)'])


%% Just look at optimal leakiness for different propA and competition values



% Responsive Strategy (Whichever Fungus is best right now) 
% Loyal Strategy (Fungus A Only) 
% Optimally Leaky (A little of both to maximize average biomass) 

% Want biomass timeseries for these three strategies across PropA values, 

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
rtot = .2;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

env_period_vals = [4*365];

results = nan(length(env_period_vals), length(leakiness_vals)); 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

propA_vals = 0:.1:1; 
competition_values = [0.01 0.04 0.05 0.1 0.2] ; 
optimal_leakiness_values = nan(length(competition_values), length(propA_vals)); 

for d = 1:length(competition_values)
    d2_2 = competition_values(d);  %density dependent mortality effect of F2 on F2
    d1_1 = d2_2;  %density dependent mortality effect of F1 on F1

for p = 1:length(propA_vals)
   propA = propA_vals(p);

   env_period = env_period_vals(j);
   envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

%Find optimal leakiness
optimal_leakinesss_sol = []; 
optimal_leakiness_biomass = 0; 
optimal_leakiness_value = []; 

leakiness_vals = 0:.05:1; 
results = nan(1, length(leakiness_vals));
for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i) = mean(final_res(1,:)); 

        if optimal_leakiness_biomass < mean(final_res(1,:))
            optimal_leakiness_biomass = mean(final_res(1,:)); 
            optimal_leakiness_sol = final_res; 
            optimal_leakiness_value = leakiness; 
        end

end

optimal_leakiness_values(d, p) = optimal_leakiness_value; 

end

end

scatter(competition_values, propA_vals, 10, optimal_leakiness_values, 'filled'); 


