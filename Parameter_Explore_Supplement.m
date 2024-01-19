%First we will show solutions for 1 fungus model? 



% set parameter values randomly 
% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns

d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 


% Set timespan and environment conditions during timespan 
tspan = [1 8000]; 

%ALWAYS ENV A 
env_period = 365; 
propA = 1; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

difference_vals = [1];
difference_val = difference_vals(4-d);
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(4) = 0; %F2
x0(5) = 1; %N

preference = 1;

rtot = 0.5; 

parameter_values = [0.001:0.005:0.01];
figure
for j = 1:length(parameter_values)

    m1 = parameter_values(j); %fungus 1 mortality
    m2 = m1; %fungus 2 mortality

    
    for i = 1:200
        if i < 100
            x0(3) = rand; %F1
        else 
          x0(3) = rand*0.001; %F1
        end


sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(preference), rtot*(preference), rtot*(1- preference), rtot*(1- preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
result = mean(final_res(1,:));

hold on 
parameter_value = parameter_values(j); 
xlabel('Mortality')
scatter(parameter_value, result, '.')

    end

end
