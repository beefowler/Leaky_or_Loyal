
%testing if results hold when variable partner quality pertains to other
%fungal traits, rather htan just uptake rate


% Constant reward strategies across environmental regimes
% Partner preference strategies vs Bet-hedging
% Mirror image fungi

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
sC = 0.005; %loss of carbon pool to environment
sN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns
u_bar = .5; %mean uptake of Nitrogen by 


u1_A = u_bar; 
u2_A = u_bar; 
u1_B = u_bar; 
u2_B = u_bar; 


%density dependent mortality parameters 
d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

%maximum reward rate
rtot = 0.8; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

% Set timespan and environment conditions during timespan 
tspan = [1 5000]; 

%environmental regimes
env_period = 365; 
propA_vals = 0:.02:1; 

%initialize results
results = nan(3,length(propA_vals)); %mean values for each strategy
results2A = nan(3, length(propA_vals)); %median and iqr for preference A
results2B = results2A; %median and iqr for preference B


figure 
clf 

    e1_A = .01;
    e1_B = .01;
    e2_A = .01;
    e2_B = .01;

for d = 1:3 %cycle through environmental severity 
    difference_vals = [0.001 0.004 0.005]; 
    difference_val = difference_vals(4-d); 
    m1 = 0.005; %fungus 1 mortality

    m1_A = m1 - difference_val;
    m1_B = m1 +difference_val;
    m2_A = m1 +difference_val;
    m2_B = m1 - difference_val;

for i = 1:length(propA_vals) %cycle through alpha/ evenness values
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) Other_Quality_Metrics_Function(t, x, g, a, s, sC, rtot, rtot, 0, 0, e1_A, e1_B, e2_A, e2_B, m1_A, m1_B, m2_A, m2_B, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 
    results2A(1,i) = median(final_res(1,:)); 
    results2A(2,i) = quantile(final_res(1,:), .25); 
    results2A(3,i) = quantile(final_res(1,:), .75); 


    if propA == .64; % Plot one simulation example for each column 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
    end


    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) Other_Quality_Metrics_Function(t, x, g, a, s, sC, 0, 0, rtot, rtot, e1_A, e1_B, e2_A, e2_B, m1_A, m1_B, m2_A, m2_B, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    results2B(1,i) = median(final_res(1,:)); 
    results2B(2,i) = quantile(final_res(1,:), .25); 
    results2B(3,i) = quantile(final_res(1,:), .75); 


    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) Other_Quality_Metrics_Function(t, x, g, a, s, sC, rtot/2, rtot/2, rtot/2, rtot/2, e1_A, e1_B, e2_A, e2_B, m1_A, m1_B, m2_A, m2_B, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

%plot median and interquartile range for prefernce strategies 
% (bet-hedging mean = median, since all constant) 
subplot(2,3,d)
plot(propA_vals, results(3,:), 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)


xlabel('Proportion of time in environment A')
ylim([0 20])
legend({'Bet-hedging'; ''; 'F_1 preference'; ''; 'F_2 preference'; ''; ''; ''})

hold on %put symbol on conditions that simulation is from 
scatter(0.64, results2A(1, propA_vals == 0.64), 25, [1 1 1], 'filled', 'MarkerEdgeColor', 'k')

end

%some formatting 
subplot(2,3,1) 
title('More severe environment')
ylabel('Photosynthetic biomass')
legend off

subplot(2,3,3)
title('More mild environment')
legend off

subplot(2,3,4) 
ylabel('Biomass')
legend off

subplot(2,3,6)
legend off

set(findall(gcf,'-property','FontSize'),'FontSize',13)