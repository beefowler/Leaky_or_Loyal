%Only want to consider end-state solutions
 
 clear all

% % set parameter values
% g = .1; %growth of plant proportional to Nitrogen pool
% a = .04; %allocation of Carbon to Carbon pool
% s = 0.01; %senesence of tree biomass
% l = 0.005; %loss of carbon pool to environment
% e1 = 0.01; %efficiency of fungus 1 carbon uptake
% e2 = 0.01; %efficiency of fungus 2 carbon uptake
% m1 = 0.004; %fungus 1 mortality
% m2 = 0.004; %fungus 2 mortality
% u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
% u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
% u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
% u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
% mN = .1; %loss of Nitrogen from tree's stores
% Ntot = 10; %total nitrogen = N + Ns

% set parameter values
g = rand(1)*.5 
a = rand(1)*.1 
s = rand(1)*.1 %senesence of tree biomass
l = rand(1)*.1 %loss of carbon pool to environment
e1 = rand(1)*.1 %efficiency of fungus 1 carbon uptake
e2 = e1; %efficiency of fungus 2 carbon uptake
m1 = rand(1)*.1 %fungus 1 mortality
m2 = m1; %fungus 2 mortality
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = rand(1)*.1 %loss of Nitrogen from tree's stores
Ntot = rand(1)*20 %total nitrogen = N + Ns


% Set timespan and environment conditions during timespan 

figure(5)
clf

    difference_val = 1;
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N
% 
d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

leakiness_vals = [0:.025:1]; 
propA_vals = [0.5:0.025:1]; 
rtot = 0.2; 

env_periods = [365/2 365 4*365] ;

%% 
for k = 1:length(env_periods); 
     env_period = env_periods(k);

     results = nan(length(leakiness_vals), length(propA_vals)); 
     coexistence_results = nan(length(leakiness_vals), length(propA_vals)); 
     for j = 1:length(propA_vals)
        
        propA = propA_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

       
for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i); 
    converged = 0; 
    tspan = [1 10000]; 

    while converged == 0

    %simulate
    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);

    thirdtolast = deval(sol, tspan(2)-env_period*3:tspan(2)-env_period*2); 
    last = deval(sol, tspan(2)-env_period:tspan(2)); 

    %check for convergence and no extinction
    coexist = 0; 
    if max(thirdtolast(1,:)) >= max(last(1,:))*.99 & min(thirdtolast(1,:)) <= min(last(1,:)*1.01) %tree biomass converging
        converged = 1; 
    else
        biomass = deval(sol, tspan(1):tspan(2));
        running_mean = movmean(biomass(1,:), 2*env_period); %if mean has changed by within 1 unit
        if range(running_mean(tspan(2)-env_period*4:tspan(2)-env_period*1))  < 1
            converged = 1;
        end
    end

    if converged == 0
        tspan(2) = tspan(2)*5
    end

    if tspan(2) > 600000
        keyboard
    end
  
    end

    if any(last(3,:)>0.01) & any(last(4,:)>0.01) %and both fungal partners are nonnegligible for some part of the cycle
            coexist = 1; 
    end

    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    biomass_val = mean(final_res(1,:)); 

    results(i,j) = biomass_val ; 
    coexistence_results(i,j) = coexist; 

    
end
end

figure(5)
subplot(2,3,k)
h = pcolor(leakiness_vals, propA_vals, results')
h.EdgeColor = 'none'; 

xlabel('Leakiness')
title(num2str(env_periods(k)))

hold on 
contour(leakiness_vals, propA_vals, coexistence_results','LevelStep', 1,'LineWidth', 2, 'color', [1 1 1], 'LineStyle', '-.')

    keyboard


end
subplot(2,3,1)
title('6 months')
ylabel('Proportion of time in Environment A')
box on 

subplot(2,3,2)
title('1 year')
box on 

subplot(2,3,3)
title('4 years')
box on 

h = colorbar; 
h.Label.String = 'Tree biomass'; 

colormap(cmocean('deep'))

%%
subplot(2,2,3)
title('Low reward rate')

colorval = cmocean('deep', 6); 
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
        plot(leakiness_vals, results, 'linewidth', 2, 'color', colorval(k+1, :));

    end
end

    subplot(2,2,3)
    title('Low reward rate')
    %ylim([10 19])
        ylabel('Tree biomass')
xlabel('Leakiness')
box on 

    legend({'1 year'; '3 years'; '5 years'; '7 years'}, 'Location', 'southwest')
    subplot(2,2,4)
        %ylim([14 23])
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
