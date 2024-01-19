% Check if there is an allee effect or if extiction on tree is unstable

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
tspan = [1 30000]; 

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
initial_vals = [0:0.00005:0.001];

rtot_vals = [.1:.05:1]; 

%initialize results
results = nan(length(rtot_vals), length(initial_vals)); 

for j = 1:length(rtot_vals)
    rtot = rtot_vals(j); 

for i = 1:length(initial_vals)
    x0(3) = initial_vals(i); %F1


sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(preference), rtot*(preference), rtot*(1- preference), rtot*(1- preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
results(j,i) = mean(final_res(1,:));
% 
% subplot(1, length(initial_vals), i)
% b = plot(sol.x, sol.y([1 3 4], :));
% b(1).Color = [89,174,159]/255;
% b(2).Color =  [196/255 118/255 165/255];
% b(3).Color = [128, 180, 232]/255;
% 
% b(1).LineWidth = 2;
% b(2).LineWidth = 2;
% b(3).LineWidth = 2;
% 
% b(2).LineStyle = '-.';
% b(3).LineStyle = '--';
% 
% legend({'Tree'; 'Fungus 1'; 'Fungus 2'})
% 
% xticks([0:10*365:30000])
% xticklabels(10*[0:30000/(365)])
% xlabel('Years')

end
end

figure
h = pcolor(rtot_vals, initial_vals, (results)');
ylabel('Initial F_1 biomass')
xlabel('Reward rate')
h.EdgeColor = 'none';
colormap('bone')


%% Do same thing but adjusting mortality? 


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
tspan = [1 30000]; 

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
initial_vals = [0:0.00005:0.001];

rtot_vals = [.1:.05:1]; 

%initialize results
results = nan(length(rtot_vals), length(initial_vals)); 

for j = 1:length(rtot_vals)
    rtot = rtot_vals(j); 

for i = 1:length(initial_vals)
    x0(3) = initial_vals(i); %F1


sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot*(preference), rtot*(preference), rtot*(1- preference), rtot*(1- preference), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
results(j,i) = mean(final_res(1,:));
% 
% subplot(1, length(initial_vals), i)
% b = plot(sol.x, sol.y([1 3 4], :));
% b(1).Color = [89,174,159]/255;
% b(2).Color =  [196/255 118/255 165/255];
% b(3).Color = [128, 180, 232]/255;
% 
% b(1).LineWidth = 2;
% b(2).LineWidth = 2;
% b(3).LineWidth = 2;
% 
% b(2).LineStyle = '-.';
% b(3).LineStyle = '--';
% 
% legend({'Tree'; 'Fungus 1'; 'Fungus 2'})
% 
% xticks([0:10*365:30000])
% xticklabels(10*[0:30000/(365)])
% xlabel('Years')

end
end

figure
h = pcolor(rtot_vals, initial_vals, (results)');
ylabel('Initial F_1 biomass')
xlabel('Reward rate')
h.EdgeColor = 'none';
colormap('bone')

