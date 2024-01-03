
% Want to compare performance of 50:50 strat with 0% leakiness

%specifically looking for parameter values where 50:50 does better

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
rtot = .4; 


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
leakiness_vals = [0 .5]; 
env_period_vals = [365/5 365/4 365/3 365/2 365 365*5];
propA_vals = 0:.1:1; 


tspan = [1 8000]; 
count = 1;
count2 = 1;

figure(1)

    difference_vals = [.5:.05:1]; 

results = nan(length(difference_vals), length(env_period_vals)); 

for d = 1:length(difference_vals)
    difference_val = difference_vals(d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;


for i = 1:length(env_period_vals)
    env_period = env_period_vals(i); 
    
        propA = 0.5; 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        leakiness = 0; 
        r1_A = (1-leakiness).*rtot;
        r1_B = leakiness*rtot; 
        r2_A = leakiness*rtot; 
        r2_B = (1-leakiness)*rtot; 

        sol1 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res1 = deval(sol1, tspan(2)-env_period*3:tspan(2));

        leakiness = 0.5; 
        r1_A = (1-leakiness).*rtot;
        r1_B = leakiness*rtot; 
        r2_A = leakiness*rtot; 
        r2_B = (1-leakiness)*rtot; 


        sol2 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res2 = deval(sol2, tspan(2)-env_period*3:tspan(2));
        results(d,i) = mean(final_res1(1,:)) - mean(final_res2(1,:));

            xticks([0:365:5000])
            xlabel('Days') 
        
    end


end

figure(5)
pcolor(results)
h = colorbar
h.Label.String = 'Cost of Maintaining Diversity'



%% (this is copy of wrapper 3 but with rtot adjustable) 

% d = 3
%     difference_vals = [.6 .8 1]; 
%     difference_val = difference_vals(4-d); 
%     u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
%     u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
%     u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
%     u2_B = difference_val;

    rtot = .2

    env_period_vals = [365/5 365/4 365/3 365/2 365 365*2 365*3 365*4 365*5 365*7];
propA_vals = 0:.1:1; 
results = nan(length(propA_vals), length(env_period_vals)); 


for i = 1:length(env_period_vals)
    env_period = env_period_vals(i); 
    
    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        leakiness = 0; 
        r1_A = (1-leakiness).*rtot;
        r1_B = leakiness*rtot; 
        r2_A = leakiness*rtot; 
        r2_B = (1-leakiness)*rtot; 

        sol1 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res1 = deval(sol1, tspan(2)-env_period*3:tspan(2));

        leakiness = 0.5; 
        r1_A = (1-leakiness).*rtot;
        r1_B = leakiness*rtot; 
        r2_A = leakiness*rtot; 
        r2_B = (1-leakiness)*rtot; 


        sol2 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res2 = deval(sol2, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res1(1,:)) - mean(final_res2(1,:));

            xticks([0:365:5000])
            xlabel('Days')
            
        
    end


end

figure
pcolor(env_period_vals, propA_vals, results)
ylabel('PropA')
xlabel('Period')
h = colorbar
h.Label.String = 'Cost of Maintaining Diversity'



%% Play with reward rates and look at simulations

% Solve ODES to run a simulation of Leaky or Loyal model 
% plot simulation results 


% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .05; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
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

rtot = 0.2; 
leakiness = 0.5; 
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

% Set timespan and environment conditions during timespan 
tspan = [1 5000]; 

env_period = 2*365; 
propA = .6; 
 envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

%run simulatio
[tout, yout] = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);


% plot

%figure
subplot(2, 1,1)
b = plot(tout, yout);
b(1).Color = [0 .9 .3]; 
b(2).Color = [0 .8 .7]; 
b(3).Color = 'r'; 
b(4).Color = 'b'; 
b(5).Color = [.5 0 .5]; 

b(2).LineStyle = '--'; 
b(5).LineStyle = '--'; 


legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
title(['rtot = ' num2str(rtot)])

rtot = 0.2; 
leakiness = 0; 
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

%run simulatio
[tout, yout] = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);

% plot

%figure
subplot(2, 1,2)
b = plot(tout, yout);
b(1).Color = [0 .9 .3]; 
b(2).Color = [0 .8 .7]; 
b(3).Color = 'r'; 
b(4).Color = 'b'; 
b(5).Color = [.5 0 .5]; 

b(2).LineStyle = '--'; 
b(5).LineStyle = '--'; 
title(['rtot = ' num2str(rtot)])


