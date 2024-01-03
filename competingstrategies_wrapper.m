% Want to have plots for different biomass for 
% 0 leakiness, 50:50, and 100% A. 



% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.05; %density dependent mortality of F1
d2 = 0.05; %density dependent mortality of F2


mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

rtot = 0.4; 


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

%initialize results
propA_vals = 0:.02:1; 
results = nan(3,length(propA_vals)); 

results2A = nan(3, length(propA_vals)); 
results2_0 = nan(3, length(propA_vals)); 


tspan = [1 8000]; 

difference_vals = [1 .9 .8 .7 .6]; 
env_period_vals = [.5 3 6]; 

for u = 1:5
    
    difference_val = difference_vals(u);
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for ep = 1:3; 
    env_period = env_period_vals(ep)*365; 


for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l,  1, 2, 0, 0, e1,  e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 
    results2A(1,i) = median(final_res(1,:)); 
    results2A(2,i) = quantile(final_res(1,:), .25); 
    results2A(3,i) = quantile(final_res(1,:), .75); 


    %next run simulation for reward strategy 0% leakiness
    %environments

    leakiness = 0; 
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 


    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l,  r1_A, r1_B, r2_A, r2_B, e1,  e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, (tspan(2)-env_period*3):tspan(2));
    results(2,i) = mean(final_res(1,:)); 
    results2_0(1,i) = median(final_res(1,:)); 
    results2_0(2,i) = quantile(final_res(1,:), .25); 
    results2_0(3,i) = quantile(final_res(1,:), .75);

    % run simulation for reward strategy 50:50

    leakiness = 0.5; 
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 


    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l,  r1_A, r1_B, r2_A, r2_B, e1,  e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, (tspan(2)-env_period*3):tspan(2));
    results(3,i) = mean(final_res(1,:)); 




end

figure(1)
subplot(3,5,(ep-1)*5+u)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
hold on 
plot(propA_vals, results2A(2,:), 'color', [0 0.4470 0.7410], 'linestyle', '--')
plot(propA_vals, results2A(3,:), 'color', [0 0.4470 0.7410], 'linestyle', '--')
plot(propA_vals, results2_0(2,:), 'color', [0.8500 0.3250 0.0980], 'linestyle', '--')
plot(propA_vals, results2_0(3,:), 'color', [0.8500 0.3250 0.0980], 'linestyle', '--')
ylim([0 23]) 




end

end


