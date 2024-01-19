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
results2A = nan(3, length(propA_vals)); 
results2B = results2A; 

env_period = 365; 

figure 
clf 

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


    if propA == .64; 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y([1 3 4], :));
           % b(1).Color = [10, 156, 0]/255;
            b(1).Color = [89,174,159]/255;
            %b(2).Color = 'k';
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;
            %b(5).Color = [214, 71, 90]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})
            %title('F_1 preference')

            xticks([0:365:3000])
            xticklabels([1:3000/365])
            xlabel('Years')
    end


    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 0, 0, 1, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 

    results2B(1,i) = median(final_res(1,:)); 
    results2B(2,i) = quantile(final_res(1,:), .25); 
    results2B(3,i) = quantile(final_res(1,:), .75); 


    %run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

subplot(2,3,d)
plot(propA_vals, results(3,:), 'k', 'linewidth', 2)
hold on 
boundedline(propA_vals, results2A(1,:), [results2A(1,:)-results2A(2,:); results2A(3,:)-results2A(1,:)]', 'linewidth', 2, 'color', [196/255 118/255 165/255], 'alpha', 'transparency', 0.5)
boundedline(propA_vals, results2B(1,:), [results2B(1,:)-results2B(2,:);results2B(3,:)- results2B(1,:)]', 'linewidth', 2, 'color', [128, 180, 232]/255, 'alpha', 'transparency', 0.5)


xlabel('Proportion of time in environment A')
%title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
ylim([0 18])
legend({'Bet hedging'; ''; 'F_1 preference'; ''; 'F_2 preference'; ''; ''})

hold on 
scatter(0.64, results2A(1, propA_vals == 0.64), 25, [1 1 1], 'filled', 'MarkerEdgeColor', 'k')

end

subplot(2,3,1) 
title('More volatile environment')
ylabel('Tree biomass')
legend off

subplot(2,3,3)
title('Less volatile environment')
legend off

subplot(2,3,4) 
ylabel('Biomass')
legend off

set(findall(gcf,'-property','FontSize'),'FontSize',13)