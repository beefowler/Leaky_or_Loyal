%Given allocatinon only in env A
%Mirror image fungi, but varying differences between them 

clear all

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
sC = 0.005; %loss of carbon pool to environment
sN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

% Set timespan and environment conditions during timespan 
tspan = [1 3000]; 

%initialize results
propA_vals = 0:.02:1; 
results = nan(1,length(propA_vals)); 
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 
fungi_resultsA = results2A;
fungi_resultsB= results2B; 
fungi_results = nan(1,length(propA_vals)); 

env_period = 365; 

figure 
clf 

u_bar = .5; %mean uptake of Nitrogen by 
u1_A = u_bar; 
u2_A = u_bar; 
u1_B = u_bar; 
u2_B = u_bar; 



    count = 0; 
    rtot = .2; 
    leakiness_vals = [0 .25 .75]; 


m1 = 0.005; %fungus 1 mortality
m1_A = m1; 
m1_B = m1; 
m2_A = m1; 
m2_B = m1; 

    e1_A = .01; %uptake of Nitrogen by fungus 1 in environment type A
    e1_B = .01; %uptake of Nitrogen by fungus 1 in environment type B
    e2_A = .01; %uptake of Nitrogen by fungus 2 in environment type A
    e2_B = .01;

 difference_val = 0.01; 
    % e1_A = .01 + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    % e1_B = .01 -difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    % e2_A = .01 -difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    % e2_B = .01 + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

    difference_val = 0.002; 
    m1_A = m1 - difference_val; 
    m1_B = m1 + difference_val; 
    m2_A = m1 + difference_val; 
    m2_B = m1 - difference_val; 


for col = 1:3
    leakiness = leakiness_vals(col)

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %next run simulation for reward strategy slightly leaky strategy 0.3 
    %environments
    sol = ode45(@(t, x) Other_Quality_Metrics_Function(t, x, g, a, s, sC, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1_A, e1_B, e2_A, e2_B, m1_A, m1_B, m2_A, m2_B, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
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


    if i == 33
         subplot(3,3,6+col)
            hold on 
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Bet hedging: Tree'; 'Nonconstant reward: Tree'; 'Nonconstant reward: Fungus 1'; 'Nonconstant reward: Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
    end
    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) Other_Quality_Metrics_Function(t, x, g, a, s, sC, rtot/2, rtot/2, rtot/2, rtot/2, e1_A, e1_B, e2_A, e2_B, m1_A, m1_B, m2_A, m2_B, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(i) = mean(final_res(1,:)); 

    fungi_results(i) = mean(final_res(3,:)) + mean(final_res(4,:)); 


    if i == 31
         subplot(3,3,6+col)
         hold on 
         b = plot(sol.x, sol.y([1], :));
         b(1).Color = 'k';
         b(1).LineWidth = 2;
         %keyboard
    end

end

subplot(3,3,col)
plot(propA_vals, results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [89,174,159]/255, 'alpha', 'transparency', 0.4)


%xlabel('Proportion of time in environment A')
%title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 20])
legend({'Bet-hedging'; ''; 'Responsive'})



subplot(3,3,col+3)
plot(propA_vals, fungi_results, 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, fungi_resultsA(1,:), [fungi_resultsA(1,:)-fungi_resultsA(2,:); fungi_resultsA(3,:)-fungi_resultsA(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, fungi_resultsB(1,:), [fungi_resultsB(1,:)-fungi_resultsB(2,:);fungi_resultsB(3,:)- fungi_resultsB(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)
xlabel('Proportion of time in environment A')
ylim([0 2])

end

subplot(3,3,1) 
title('No leakiness (0)')
ylabel('Tree biomass')

subplot(3,3,2)
title('Moderate leakiness (.25)')

subplot(3,3,3)
title('High leakiness (.75)')

subplot(3,3,4) 
ylabel('Fungal biomass')
legend({'Bet-hedging: both fungi'; ''; 'Responsive: Fungus 1'; ''; 'Responsive: Fungus 2'})


set(findall(gcf,'-property','FontSize'),'FontSize',13)