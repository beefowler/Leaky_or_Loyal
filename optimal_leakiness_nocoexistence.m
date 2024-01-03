

envPeriod_vals = 356*[3] ;
propA_vals = [0:.02:1];

optimal_leakiness_values = nan(length(envPeriod_vals), length(propA_vals));
optimal_leakiness_biomass = optimal_leakiness_values;


for p = 1:length(envPeriod_vals)

    for a = 1:length(propA_vals)

        [X2,FVAL2,~,~] = fmincon(@(x) negbiomass(x, envPeriod_vals(p), propA_vals(a)), 1, [],[],[],[],0,1);

        %
        % Xvec = [X1 X2 X3];
        % Fvec = [FVAL1 FVAL2 FVAL3];
        % [~,a] = min(Fvec);
        % optimal_leakiness_values(d, 1) = Xvec(a);
        % optimal_leakiness_values(d, 2) = Fvec(a);

        optimal_leakiness_values(p, a) = X2;
        optimal_leakiness_biomass(p, a) = -FVAL2;


    end

    figure(4)
    subplot(1,2,1)
    hold on 
    plot(propA_vals, optimal_leakiness_values(p,:))
    ylabel('Optimal leakiness')
    subplot(1,2,2)
    hold on 
    plot(propA_vals, optimal_leakiness_biomass(p,:))

    ylabel('Optimal mean biomass')
    xlabel('Proportion of time in environment A')


end

legend({'1/2 year'; '2 years'; '8 years'})



function B = negbiomass(leakiness, env_period, propA)
% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns
rtot = 0.2;

difference_val = 1; 
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;

r1_A = (1-leakiness).*rtot;
r1_B = leakiness*rtot;
r2_A = leakiness*rtot;
r2_B = (1-leakiness)*rtot;

x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


tend = max(3000, env_period.*6);
tspan = [1 tend];

envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

B = -mean(final_res(1,:)); %negative mean biomass


end

