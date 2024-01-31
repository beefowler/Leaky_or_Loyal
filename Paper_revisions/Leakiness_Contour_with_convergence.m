
% Code to generate Figure 3 in Stevens et al.

% Want biomass for varying alpha /evenness and leakiness.
% Distinct panels for different periods.
% Only want to consider end-state solutions, so make sure solutions
% converge, even for long period cycles.

clear all
figure(5)
clf

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
sC = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.005; %fungus 1 mortality
m2 = 0.005; %fungus 2 mortality
sN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

u_bar = .5; %mean uptake of nitrogen by fungus
difference_val = .5; %severity of environment, high
u1_A = u_bar + difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = u_bar - difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = u_bar - difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = u_bar + difference_val;%uptake of Nitrogen by fungus 2 in environment type B

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

%density dependent mortality paramters 
d1_1 = .1;
d2_1 = .05;
d2_2 = .1;
d1_2 = .05;

leakiness_vals = [0:.025:1];
propA_vals = [0.5:0.025:1];
rtot = 0.2;

%three panels, three environmental periods
env_periods = [365/2 365 5*365] ;

%%
for k = 1:length(env_periods) % For each of the panels
    env_period = env_periods(k);

    %initialize results for this period 
    results = nan(length(leakiness_vals), length(propA_vals));
    coexistence_results = nan(length(leakiness_vals), length(propA_vals)); %check if both fungus still exist

    for j = 1:length(propA_vals) % vary evenness

        propA = propA_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;


        for i = 1:length(leakiness_vals) % vary leakiness values
            leakiness = leakiness_vals(i);
            converged = 0;
            tspan = [1 10000];

            while converged == 0 %keep repeating until solution has converged

                %simulate
                sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, sC, rtot*(1-leakiness), rtot*leakiness, rtot*leakiness, rtot*(1-leakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, sN, Ntot, envA_treat), tspan, x0);

                thirdtolast = deval(sol, tspan(2)-env_period*3:tspan(2)-env_period*2);
                last = deval(sol, tspan(2)-env_period:tspan(2));

                %check for convergence and no extinction
                coexist = 0;
                if max(thirdtolast(1,:)) >= max(last(1,:))*.99 && min(thirdtolast(1,:)) <= min(last(1,:)*1.01) %tree biomass converging
                    converged = 1;
                else
                    biomass = deval(sol, tspan(1):tspan(2));
                    running_mean = movmean(biomass(1,:), 2*env_period); %if mean has changed by within 1 unit
                    if range(running_mean(tspan(2)-env_period*4:tspan(2)-env_period*1))  < 1
                        converged = 1;
                    end
                end

                if converged == 0
                    tspan(2) = tspan(2)*5 %spit out timespan so I know how hard it's working
                end

                if tspan(2) > 600000
                    keyboard %stop if goes on too long, check what's goingo n
                end

            end

            if any(last(3,:)>0.01) && any(last(4,:)>0.01) %and both fungal partners are nonnegligible for some part of the cycle
                coexist = 1;
            end

            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:)); %final biomass is from last three periods

            results(i,j) = biomass_val ;
            coexistence_results(i,j) = coexist;


        end
    end

    figure(5)
    subplot(1,3,k)
    h = pcolor(leakiness_vals, propA_vals, results');
    h.EdgeColor = 'none';

    xlabel('Leakiness')
    title(num2str(env_periods(k)))

    hold on
    contour(leakiness_vals, propA_vals, coexistence_results','LevelStep', 1,'LineWidth', 2, 'color', [1 1 1], 'LineStyle', '-.')

end

%some more formating and titles
subplot(1,3,1)
title('T = 6 months')
ylabel('Proportion of time in Environment A')
box on

subplot(1,3,2)
title('T = 1 year')
box on

subplot(1,3,3)
title('T = 5 years')
box on

h = colorbar;
h.Label.String = 'Tree biomass';

set(findall(gcf,'-property','FontSize'),'FontSize',13)

%using colormap from https://matplotlib.org/cmocean/
colormap(cmocean('tempo'))
