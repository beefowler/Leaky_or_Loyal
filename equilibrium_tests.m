
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
rtot = .2;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

env_period= 1; 
propA = 1;
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

r1_A = rtot;
r1_B = rtot;
r2_A = rtot;
r2_B = rtot;

uval = 0.6;
u1_A = uval;
u1_B = uval;
u2_A = uval;
u2_B = uval;


sol1 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);

figure
b = plot(sol1.x, sol1.y);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';

%% Now let's try with a coexistence possible equilibriumm

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.2;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.2;  %density dependent mortality effect of F2 on F2
rtot = .2;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns



% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

env_period= 1; 
propA = 1;
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;

rdiff = .1; 
r1_A = rtot-rdiff;
r1_B = rtot-rdiff;
r2_A = rtot;
r2_B = rtot;

difference_val = .5;
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;

sol1 = ode45(@(t, x)  leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);



figure
b = plot(sol1.x, sol1.y);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';


%% Nice those parameters allow for coexistence at relative levels that are modifiable by e and r 


difference_val = 1;
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;


%initialize results
propA_vals = 0:.1    :1; 
results = nan(3,length(propA_vals)); 
env_period = 365; 

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 66% fungus 1 in both environments
    sol = ode45(@(t, x)  leaky_or_loyal_coexistence(t, x, g, a, s, l, .2, .2, .1, .1, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 


       b = plot(sol.x, sol.y);
    b(1).Color = [0 .9 .3];
    b(2).Color = [0 .8 .7];
    b(3).Color = 'r';
    b(4).Color = 'b';
    b(5).Color = [.5 0 .5];

    b(2).LineStyle = '--';
    b(5).LineStyle = '--';
    pause 


    %next run simulation for reward strategy 100% fungus 1 in both
    %environments
    sol = ode45(@(t, x)  leaky_or_loyal_coexistence(t, x, g, a, s, l, .3, .3, 0, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 


    % run simulation for reward strategy 50:50
    sol = ode45(@(t, x)  leaky_or_loyal_coexistence(t, x, g, a, s, l, .15, .15, .15, .15, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

end

figure
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
legend({'66%A'; '100%A'; '50%A'})


%%

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.2;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.2;  %density dependent mortality effect of F2 on F2
rtot = .8;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


leakiness_vals = 0:.1:1; 
env_period_vals = [180 365*2 365*4];
propA_vals = [0.1 .4 .7  1]; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

for p = 1:4
    propA = propA_vals(p); 

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

    for j = 1:length(env_period_vals)

        env_period = env_period_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));
    end


end

figure(2)
subplot(1,4,p)
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend({'6 months'; '2 years'; '4 years'}); 
title(['PropA = ' num2str(propA)])

end



%% 
leakiness_vals = 0:.05:1; 
env_period =  365*2;
propA = 0.4; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

plot_leakiness = [0 0.1 leakiness_vals(8) 0.5 0.8 1];

envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 
count = 1; 

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

    sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(j,i) = mean(final_res(1,:));
 
    figure(1)
    if ~isempty(intersect(leakiness, plot_leakiness)) 
            subplot(2,6,3+count)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
             xlabel('Days')
             title(['leakiness = ' num2str(leakiness)])
            count = count+1;
            if count == 4
                count = 7 
            end

    end
end

figure(1)
subplot(1,2,1)
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend({'6 months'; '2 years'; '4 years'}); 
title(['PropA = ' num2str(propA)])


%% Test original model for sensitivity to initial conditions :( 


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


%initialize results
leakiness_vals = 0:.05:1; 
env_period_vals = [180 365 365*2 365*3];
propA = .5; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

tspan = [1 8000]; 

rtot = .2; 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 10; %F1
x0(4) = 2; %F2
x0(5) = 1; %N


for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

    for j = 1:length(env_period_vals)

        env_period = env_period_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

    end


    if leakiness == 0.6;
        figure
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
             xlabel('Days')
             title(['leakiness = ' num2str(leakiness)])
            count = count+1;
            pause
    end

end

figure
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend({'6 months'; '1 year'; '2 years'; '3 years'}); 

%% 
% If my understanding is correct, new version should not be sensitive to initial conditons

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1_1 = 0.2;  %density dependent mortality effect of F1 on F1
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d2_2 = 0.2;  %density dependent mortality effect of F2 on F2
rtot = .8;
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


leakiness_vals = 0:.05:1; 
env_period_vals = [2*365];
propA = .8; 
results = nan(length(env_period_vals), length(leakiness_vals)); 

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

tspan = [1 8000];

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

    for j = 1:length(env_period_vals)

        env_period = env_period_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        
        %sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));
    end

end

figure
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
title(['PropA = ' num2str(propA)])




