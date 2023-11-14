%% 
% We want to now test whether the 50:50 leakiness is still optimal if the
% tree has the ability to change it's allocation amount. H: Probably not. 



% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a_A = .05; %allocation of Carbon to Carbon pool in environment type A
a_B = .05; %allocation of Carbon to Carbon pool in environment type B
s = 0.001; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
r1_A = 0; %reward rate to fungus 1 in environment type A
r1_B = 0.2; %reward rate to fungus 1 in environment type B
r2_A = 0.2; %reward rate to fungus 2 in environment type A
r2_B = 0; %reward rate to fungus 2 in environment type B
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
u1_A = .6; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = .4; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = .4; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = .6; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
propA_vals = 0:.02:1; 
results = nan(3,length(propA_vals)); 

tspan = [1 20000]; 
env_period = 1000;


for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 


    if propA == 0.04; 
          figure
            subplot(2, 1, 1)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})

            subplot(2,1,2)
            plot(sol.x, envA_treat(sol.x))
            yticks([0 1])
            yticklabels({'B'; 'A'})
            title('Environment A')
    end


    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B,  s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 


    %next run simulation for reward strategy 50:50
    %environments
    sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 


end

figure
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')

%%

%Mirror image fungi, but varying differences between them 

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a_A = .05; %allocation of Carbon to Carbon pool in environment type A
a_B = .05; %allocation of Carbon to Carbon pool in environment type B
s = 0.001; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
r1_A = 0; %reward rate to fungus 1 in environment type A
r1_B = 0.2; %reward rate to fungus 1 in environment type B
r2_A = 0.2; %reward rate to fungus 2 in environment type A
r2_B = 0; %reward rate to fungus 2 in environment type B
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


u1_A = .5; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = .5; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = .5; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = .5; %uptake of Nitrogen by fungus 1 in environment type B

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 

tspan = [1 4000]; 
env_period = 365; 

for a = 1:6; 
    allocation_diffs = [.5:.1:1];
    a_A = allocation_diffs(7-a).*0.05; 
    a_B = (1-allocation_diffs(7-a)).*0.05; 

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 



    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot,  envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    % % run simulation for reward strategy 50:50
    % sol = ode45(@(t, x) leaky_or_loyal_variable_a(t, x, g, a_A, a_B, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot,  envA_treat), tspan, x0);
    % final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    % results(3,i) = mean(final_res(1,:)); 

if i == 23

figure(1)
subplot(2,3,a)
% plot(propA_vals, results)
% ylabel('Average tree biomass')
% xlabel('Proportion of time in environment A')
%legend({'100%A'; '100%B'; '50:50'})

 b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
title(['difference: ' num2str(allocation_diffs(7-a))])
end



end


figure(2)
subplot(2,3,a)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')


end


