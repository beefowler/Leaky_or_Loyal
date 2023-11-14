% Solve ODES to run a simulation of Leaky or Loyal model 
% plot simulation results 


% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .05; %allocation of Carbon to Carbon pool
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

envA = @(t) mod(floor(t/(365/2)), 2); %365 day period with half A and half B

%run simulatio
[tout, yout] = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA), tspan, x0);


% plot

figure
subplot(2, 1,1)
b = plot(tout, yout);
b(1).Color = [0 .9 .3]; 
b(2).Color = [0 .8 .7]; 
b(3).Color = 'r'; 
b(4).Color = 'b'; 
b(5).Color = [.5 0 .5]; 

b(2).LineStyle = '--'; 
b(5).LineStyle = '--'; 


legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
title('Simulation')


subplot(2,1,2)
plot(tout, envA(tout))
yticks([0 1])
yticklabels({'B'; 'A'})
title('Environment A')

%Alright this reproduces Holly's results from the R script. Looking good. 

%% Let's start running experiments

%First we want to test how biomass responds to changes in enviornmental
%variability for 3 different reward strategies. 

%need a better way to adjust environment function: 

env_period = 365; 
prop_A = 2/3; 
%this should divide timeline into periodic chunks with first prop_A portion
%A and then latter part of period portion B (0)
envA_treat = @(t) discretize(rem(t, env_period), [0 prop_A*env_period env_period]) == 1 ; 

%but remember you need to redefine envA_treat every time you change prop_A
%and env_period variable values


%% 

%initialize results
propA_vals = 0:.02:1; 
results = nan(3,length(propA_vals)); 

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

    % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

figure
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')


%% now vary reward strategy in fixed environmental regime

%initialize results
leakiness_vals = 0:.1:1; 
env_period_vals = [180 365 365*2 365*3];
propA = .5; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

tspan = [1 8000]; 

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

figure
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend({'6 months'; '1 year'; '2 years'; '3 years'}); 


%%

%Same thing but with un-equal fungi


u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B

%initialize results
leakiness_vals = 0:.02:1; 
propA = .5; 

%uptake values for B
u2_Bvals = 0:.1:1;
results = nan(length(u2_Bvals), length(leakiness_vals)); 

tspan = [1 8000]; 

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1-leakiness; 

    for j = 1:length(u2_Bvals)

        u2_B = u2_Bvals(j); 

        env_period = 365;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));
    end

end

figure
colormat = viridis(11);
for j = 1:length(u2_Bvals)
    plot(leakiness_vals, results(j,:), 'color', colormat(j,:))
    hold on 
end
ylabel('Average tree biomass')
xlabel('Leakiness')


%% 

%Mirror image fungi, but varying differences between them 

u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 

env_period = 365; 

for d = 1:6; 
    difference_vals = .5:.1:1
    difference_val = difference_vals(7-d); 
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

    % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

subplot(2,3,d)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
title(['difference: ' num2str(difference_val)])
ylim([0 18])
legend({'100%A'; '100%B'; '50:50'})
end


%%

% We now need to vary leakiness, proportion of time in environment A


%Mirror image Fungi again
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = .6; %uptake of Nitrogen by fungus 1 in environment type B

%initialize results
propA_vals = 0:.02:1; 
leakiness_vals = 0:.05:.5; 
results = nan(length(leakiness_vals), length(propA_vals)); 

tspan = [1 8000]; 
env_period = 800; 

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    for j = 1:length(leakiness_vals)
        leakiness = leakiness_vals(j);

        r1_A = 1-leakiness;
        r1_B = leakiness;
        r2_A = leakiness;
        r2_B = 1-leakiness;

        %first run simulation for reward strategy 100% fungus 1 in both environments
        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

        if propA == .6
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
            title(num2str(leakiness))


            subplot(2,1,2)
            plot(sol.x, envA_treat(sol.x))
            yticks([0 1])
            yticklabels({'B'; 'A'})
            title('Environment A')
        end

    end

end

figure
colormat = viridis(length(leakiness_vals));
for j = 1:length(u2_Bvals)
    plot(propA_vals, results(j,:), 'color', colormat(j,:))
    hold on 
end
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')









