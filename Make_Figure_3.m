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
results = nan(1,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 
fungi_resultsA = results2A;
fungi_resultsB= results2B; 
fungi_results = nan(1,length(propA_vals)); 

env_period = 365; 

figure 
clf 

    difference_val = 1;
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;


    count = 0; 
    rtot = .2; 
    leakiness_vals = [0 .25 .75]; 

for col = 1:3
    leakiness = leakiness_vals(col)

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    % %first run simulation for reward strategy 0% leaky
    % sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    % final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    % results(1,i) = mean(final_res(1,:)); 
    % results2A(1,i) = median(final_res(1,:)); 
    % results2A(2,i) = quantile(final_res(1,:), .25); 
    % results2A(3,i) = quantile(final_res(1,:), .75); 
    % 


    %next run simulation for reward strategy slightly leaky strategy 0.3 
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

    results2B(1,i) = median(final_res(1,:)); 
    results2B(2,i) = quantile(final_res(1,:), .25); 
    results2B(3,i) = quantile(final_res(1,:), .75); 

     fungi_resultsA(1,i) = median(final_res(3,:)); 
    fungi_resultsA(2,i) = quantile(final_res(3,:), .25); 
    fungi_resultsA(3,i) = quantile(final_res(3,:), .75); 

     fungi_resultsB(1,i) = median(final_res(4,:)); 
    fungi_resultsB(2,i) = quantile(final_res(4,:), .25); 
    fungi_resultsB(3,i) = quantile(final_res(4,:), .75); 


    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(i) = mean(final_res(1,:)); 

    fungi_results(i) = mean(final_res(3,:)); 



end

subplot(2,3,col)
plot(propA_vals, results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [64,131,159]/255, 'alpha', 'transparency', 0.4)


%xlabel('Proportion of time in environment A')
%title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 18])
legend({'Bet hedging'; ''; 'Switching strategy'})



subplot(2,3,col+3)
plot(propA_vals, fungi_results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, fungi_resultsA(1,:), [fungi_resultsA(1,:)-fungi_resultsA(2,:); fungi_resultsA(3,:)-fungi_resultsA(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, fungi_resultsB(1,:), [fungi_resultsB(1,:)-fungi_resultsB(2,:);fungi_resultsB(3,:)- fungi_resultsB(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
xlabel('Proportion of time in environment A')
ylim([0 .8])

end

subplot(2,3,1) 
title('No leakiness (0)')
ylabel('Tree biomass')

subplot(2,3,2)
title('Moderate leakiness (.25)')

subplot(2,3,3)
title('High leakiness (.75)')

subplot(2,3,4) 
ylabel('Fungal biomass')
legend({'Bet hedging'; ''; 'Switching strategy: Fungus 1'; ''; 'Switching strategy: Fungus 2'})


set(findall(gcf,'-property','FontSize'),'FontSize',13)