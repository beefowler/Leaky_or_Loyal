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
x0(3) = 3; %F1
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


leakiness_vals = [0:.05:1]; 
propA_vals = [0.5:0.1:1]; 
rtot = 0.2; 

env_period = 365;

results = nan(length(propA_vals), length(leakiness_vals)); 

comp_vals = [0.05 0.075 0.1 0.125 0.15];
colorvals = viridis(length(comp_vals)); 
for k = 1:length(comp_vals)
    d1_1 = .1; 
    d2_1 = comp_vals(k); 
    d2_2 = .1; 
    d1_2 = comp_vals(k); 

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

    results(j,i) = biomass_val ; 
   
end

subplot(1,length(propA_vals), j)
hold on
plot(leakiness_vals, results(j, :), 'linewidth', 2, 'color', colorvals(k,:))
ylabel('Tree biomass')
xlabel('leakiness')

end





end








