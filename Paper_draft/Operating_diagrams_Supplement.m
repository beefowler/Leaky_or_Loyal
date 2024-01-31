%Parameter space explore supplement

%keep environment constant 

%Vary tree growth rate and tree sensence 
%Vary density independent mortality and carbon conversion efficeincy 
%Vary carbon convesrion efficiency and allocation rate 
%vary total nitrogen adn rate of nitrogen loss from tree 


% Initial conditions and parameters we wont change
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

rtot = 0.5;

d1_1 = .1; 
d2_1 = .05; 
d2_2 = .1; 
d1_2 = .05; 

u1_A = 1; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 0; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 0; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = 1; %uptake of Nitrogen by fungus 1 in environment type B

tspan = [1 3000]; 

env_period = 365; 
propA = .5; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 


%% Nitrogen and rate of Nitrogen loss

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality


parameter1_vals = [1:5:50];
parameter2_vals = [0.001 0.01 0.1 1 10];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    Ntot = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       mN = parameter2_vals(j); 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 
    end

end

figure

subplot(1,2,1)
imagesc([1:length(parameter2_vals)], 'Ydata', parameter1_vals, 'Cdata', results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
ylabel({'N_T + N_S'})
xlabel('m_N')
set(gca, 'ydir', 'normal')


subplot(1,2,2)
imagesc([1:length(parameter2_vals)], 'Ydata', parameter1_vals, 'Cdata', var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
xlabel('m_N')
set(gca, 'ydir', 'normal')


colormap(cmocean('tempo'))

%% growth of plant and senescence of tree 

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
mN = .1; 
Ntot = 10; 

parameter1_vals = [0.001 0.01 0.1 1 10];
parameter2_vals = [0.001 0.01 0.1 1 10];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    g = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       s = parameter2_vals(j); 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 
    end

end

figure

subplot(1,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
ylabel('g')
xlabel('s')
set(gca, 'ColorScale', 'log')
set(gca, 'ydir', 'normal')


subplot(1,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('s')
set(gca, 'ydir', 'normal')


colormap(cmocean('tempo'))
set(gca, 'ColorScale', 'log')


%% growth of plant and nitrogen loss

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
mN = .1; 
Ntot = 10; 

parameter1_vals = [0.001 0.01 0.1 1 10];
parameter2_vals = [0.001 0.01 0.1 1 10];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    g = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       mN = parameter2_vals(j); 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 
    end

end

figure

subplot(1,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
ylabel('g')
xlabel('m_N')
set(gca, 'ColorScale', 'log')
set(gca, 'ydir', 'normal')


subplot(1,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m_N')
set(gca, 'ydir', 'normal')


colormap(cmocean('tempo'))
set(gca, 'ColorScale', 'log')


%% carbon allocation rate and carbon conversion efficiency

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
mN = .1; 
Ntot = 10; 

parameter1_vals = [0.001 0.005 0.01 0.05 0.1];
parameter2_vals = [0.001 0.005 0.01 0.05 0.1];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    a = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       e1 = parameter2_vals(j); 
        e2 = e1;

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 
    end

end

figure

subplot(1,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
ylabel('a')
xlabel('e1 & e2')
set(gca, 'ydir', 'normal')


subplot(1,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('e1 & e2')
set(gca, 'ydir', 'normal')



colormap(cmocean('tempo'))

%% carbon conversion efficiency and fungal mortality

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.004; %fungus 1 mortality
m2 = 0.004; %fungus 2 mortality
mN = .1; 
Ntot = 10; 

parameter1_vals = [0.001 0.005 0.01 0.05];
parameter2_vals = [0.001 0.005 0.01 0.05];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    e1 = parameter1_vals(i);
    e2 = e1; 

    for j = 1:length(parameter2_vals)
       m1 = parameter2_vals(j); 
        m2 = m1;

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 
    end

end

figure

subplot(1,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m1 & m2')
ylabel('e1 & e2')
set(gca, 'ydir', 'normal')


subplot(1,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m1 & m2')
set(gca, 'ydir', 'normal')


colormap(cmocean('tempo'))


