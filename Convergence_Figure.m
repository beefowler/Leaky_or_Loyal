%Only want to consider end-state solutions
 
 clear all

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.005; %fungus 1 mortality
m2 = 0.005; %fungus 2 mortality
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
tspan = [1 8000]; 

figure(5)
clf

    difference_val = 1;
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;



dvals = [0.05 0.1 0.15]; 
d1_1 = .05; 
d2_1 = .01; 
d2_2 = .05; 
d1_2 = .01; 

leakiness_vals = [0:.05:1]; 
propA_vals = [0:0.05:1]; 
rtot = 0.2; 

env_periods = [365/2 365 4*365] ;

for k = 1:length(env_periods); 
     env_period = env_periods(k);

     results = nan(length(leakiness_vals), length(propA_vals)); 
     for j = 1:length(propA_vals)
        
        propA = propA_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

       
for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i); 

       
    %simulate
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);

    thirdtolast = deval(sol, tspan(2)-env_period*3:tspan(2)-env_period*2); 
    last = deval(sol, tspan(2)-env_period:tspan(2)); 

    %check for convergence and no extinction
    converged = 0; 
    coexist = 0; 
    if max(thirdtolast(1,:)) >= max(last(1,:))*.99 & min(thirdtolast(1,:)) <= min(last(1,:)*1.01) %tree biomass converging
        converged = 1; 
    end
    if any(last(3,:)>0.01) & any(last(4,:)>0.01) %and both fungal partners are nonnegligible for some part of the cycle
            coexist = 1; 
    end

    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    biomass_val = mean(final_res(1,:)); 

    results(i,j) = biomass_val ; 

    subplot(2,length(env_periods),k)
    hold on 
    if converged & coexist 
        scatter(leakiness, propA, 40, biomass_val, 'filled', 'MarkerEdgeColor', 'k')
    elseif converged & ~coexist
        scatter(leakiness, propA, 20, biomass_val, 'd', 'filled', 'MarkerEdgecolor', 'k')
    elseif ~converged & ~coexist 
        scatter(leakiness, propA, 20, biomass_val, 'd')
    elseif ~converged & coexist 
         scatter(leakiness, propA, 40, biomass_val)
    end

   
end

% 
% if k == 1 | k == 3
%     subplot(2,2,2+k)
%     hold on 
%     if j == 11 
%         plot(leakiness_vals, results(:,j), 'linewidth', 2, 'color', 'k'); 
%     elseif j == 16
%         plot(leakiness_vals, results(:,j), '--', 'linewidth', 2, 'color','k'); 
%     elseif j == 21
%         plot(leakiness_vals, results(:,j), '-.', 'linewidth', 2, 'color', 'k'); 
%     end
%     ylabel('Tree biomass')
%     xlabel('Leakiness')
% 
% end

     end


subplot(2,3,k)
xlabel('Leakiness')
title(num2str(env_periods(k)))
caxis([0 25]); 

end
subplot(2,3,1)
title('6 months')
ylabel('Proportion of time in Environment A')

subplot(2,3,2)
title('1 year')

subplot(2,3,3)
title('4 years')

h = colorbar; 
h.Label.String = 'Tree biomass'; 

colormap viridis

%%
subplot(2,2,3)
title('Low reward rate')

colorval = viridis(4); 
rtot_vals = [.2 .8]; 
env_periods = [365 365*3 365*5 365*7];
for rt = 1:2;
    rtot = rtot_vals(rt) ;

    for k = 1:length(env_periods);
        env_period = env_periods(k);
        propA = .5;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;
        
        results = nan(1, length(leakiness_vals));
        for i = 1:length(leakiness_vals)
            leakiness = leakiness_vals(i);

            tspan = [1 8000];
            if tspan(2)-env_period*4 < 8000
                tspan = [1 env_period*8];
            end

            %simulate
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:));

            results(i) = biomass_val ;

        end
        subplot(2,2,2+rt)
        hold on
        plot(leakiness_vals, results, 'linewidth', 2, 'color', colorval(length(env_periods)+1-k, :));

    end
end

    subplot(2,2,3)
    title('Low reward rate')
    ylim([16 23])
        ylabel('Tree biomass')
xlabel('Leakiness')
box on 
    legend({'1 year'; '3 years'; '5 years'; '7 years'}, 'Location', 'southwest')
    subplot(2,2,4)
        ylim([16 23])
    xlabel('Leakiness')
    title('High reward rate')
box on 

    %
    % %
    % figure(4)
    % b = plot(sol.x, sol.y([1 3 4], :));
    % b(1).Color = [10, 156, 0]/255;
     % %b(2).Color = 'k';
     % b(2).Color = [196/255 118/255 165/255];
     % b(3).Color = [128, 180, 232]/255;
     % %b(5).Color = [214, 71, 90]/255;
     % 
     % b(1).LineWidth = 2;
     % b(2).LineWidth = 2;
     % b(3).LineWidth = 2;
     % 
     % b(2).LineStyle = '-.';
     % b(3).LineStyle = '--';
     % 
     % 
