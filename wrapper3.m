%After meeting with Holly, let's answer some questions

%Can we get fungi to stay lower than their carrying capacity for longer?

%Given allocatinon only in env A
%Mirror image fungi


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

leakiness = .2;
    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1- leakiness; 

%baseline uptake rate 
u = 1; 

mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 3000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 

env_period = 365; 

figure

for d = 1:3 
    difference_vals = [.8 .9 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = (difference_val).*u; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = (1-difference_val).*u; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = (1-difference_val).*u; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = (difference_val).*u;

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 1, 0, 0, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 

    if propA == .5; 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
           % title(['Leakiness: 0' num2str(leakiness) ', PropA: 0.8'])
        title(['Leakiness: 0' ', PropA: 0.5'])

            xticks([0:365:8000])
            xlabel('Days')
    end
    % 
    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 
        
    % if propA == .8; 
    %     subplot(2,3,3+d)
    %         b = plot(sol.x, sol.y);
    %         b(1).Color = [0 .9 .3];
    %         b(2).Color = [0 .8 .7];
    %         b(3).Color = 'r';
    %         b(4).Color = 'b';
    %         b(5).Color = [.5 0 .5];
    % 
    %         b(2).LineStyle = '--';
    %         b(5).LineStyle = '--';
    % 
    %         legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
    %         title(['Leakiness: ' num2str(leakiness) ', PropA: 0.8'])
    % 
    %         xticks([0:365:8000])
    %         xlabel('Days')
    % end

  
     % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

     

end

subplot(2,3,d)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
%ylim([0 24])
legend({'0 Leakiness'; 'Leaky'; '50:50'})
end

%% Test that that result holds for other parameters
% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1 = 0.05; %density dependent mortality of F1
d2 = 0.05; %density dependent mortality of F2

leakiness = .2;
    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1- leakiness; 

%baseline uptake rate 
u = 1; 

mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


% Set timespan and environment conditions during timespan 
tspan = [1 8000]; 

%initialize results
propA_vals = 0:.02:1; 
u2_Bvals = 0:.1:1;
results = nan(3,length(propA_vals)); 

env_period = 365; 

figure

for d = 1:3 
    difference_vals = [.8 .9 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = (difference_val).*u; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = (1-difference_val).*u; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = (1-difference_val).*u; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = (difference_val).*u;

for i = 1:length(propA_vals)
    propA = propA_vals(i); 
    envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    %first run simulation for reward strategy 100% fungus 1 in both environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 1, 0, 0, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(1,i) = mean(final_res(1,:)); 

    if propA == .8 
        subplot(2,3,3+d)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
           % title(['Leakiness: 0' num2str(leakiness) ', PropA: 0.8'])
        title(['Leakiness: ' num2str(leakiness) ', PropA:' num2str(propA)])

            xticks(0:365:8000)
            xlabel('Days')
    end
    
    %next run simulation for reward strategy 100% fungus 2 in both
    %environments
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B,  e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(2,i) = mean(final_res(1,:)); 
        
    % if propA == .8; 
    %     subplot(2,3,3+d)
    %         b = plot(sol.x, sol.y);
    %         b(1).Color = [0 .9 .3];
    %         b(2).Color = [0 .8 .7];
    %         b(3).Color = 'r';
    %         b(4).Color = 'b';
    %         b(5).Color = [.5 0 .5];
    % 
    %         b(2).LineStyle = '--';
    %         b(5).LineStyle = '--';
    % 
    %         legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
    %         title(['Leakiness: ' num2str(leakiness) ', PropA: 0.8'])
    % 
    %         xticks([0:365:8000])
    %         xlabel('Days')
    % end

  
     % run simulation for reward strategy 50:50
    sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
    final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
    results(3,i) = mean(final_res(1,:)); 

     

end

subplot(2,3,d)
plot(propA_vals, results)
ylabel('Average tree biomass')
xlabel('Proportion of time in environment A')
title(['Uptake rates: ' num2str(u1_A) ' and ' num2str(u2_A)])
%ylim([0 24])
legend({'0 Leakiness'; 'Leaky'; '50:50'})
end
%%

% Flip so axis is leakinesss and colors are propA

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1 = 0.05; %density dependent mortality of F1
d2 = 0.05; %density dependent mortality of F2

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
leakiness_vals = 0:.05:1; 
env_period = 365*2;
propA_vals = 0:.1:1; 

results = nan(length(propA_vals), length(leakiness_vals)); 

tspan = [1 8000]; 
count = 1;
count2 = 1;

figure(1)

d = 1;
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1-leakiness; 

    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

         if mod(j,2) == 1 && leakiness == 0.2 && j < 6
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
            title(['propA = ' num2str(propA) '; leakiness = 0.2'])

            xticks(0:365:8000)
            xlabel('Days')
            count = count+1;
            
         end
         if mod(j,2) == 1 && leakiness == 0.8 && j < 6
            subplot(2,6,9+count2)
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
            title(['propA = ' num2str(propA) '; leakiness = 0.8'])

            xticks(0:365:8000)
            xlabel('Days')
            count2 = count2+1;

         end

    end


end

figure(1)
subplot(1,2,1)
plot(leakiness_vals, results)
ylabel('Average tree biomass')
xlabel('Leakiness')
legend(string(propA_vals))

%% 
%let's look at simulations across one line


% Flip so axis is leakinesss and colors are propA

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d1 = 0.05; %density dependent mortality of F1
d2 = 0.05; %density dependent mortality of F2

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
leakiness_vals = 0:.05:1; 
env_period = 365*2;
propA_vals = 0:.1:1; 

results = nan(length(propA_vals), length(leakiness_vals)); 

tspan = [1 8000]; 
count = 1;
count2 = 1;

figure(1)

d = 1;
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = 1-leakiness;
    r1_B = leakiness; 
    r2_A = leakiness; 
    r2_B = 1-leakiness; 

    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

         if propA == 0.9
             if leakiness == 0
             subplot(2,6,4)
             elseif leakiness == leakiness_vals(7)
             subplot(2,6,5)
             elseif leakiness == 0.5
             subplot(2,6,6)
             elseif leakiness == 0.55
             subplot(2,6,10)
             elseif leakiness == 0.8
             subplot(2,6,11)
             elseif leakiness == 1 
             subplot(2,6,12)
             else
                 continue
             end
            b = plot(sol.x, sol.y);
            b(1).Color = [0 .9 .3];
            b(2).Color = [0 .8 .7];
            b(3).Color = 'r';
            b(4).Color = 'b';
            b(5).Color = [.5 0 .5];

            b(2).LineStyle = '--';
            b(5).LineStyle = '--';

            legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
            title(['leakiness = ' num2str(leakiness)])

            xticks([0:365:8000])
            xlabel('Days')
            
         end
    end


end

figure(1)
subplot(1,2,1)
plot(leakiness_vals, results(1:6, :))
ylabel('Average tree biomass')
xlabel('Leakiness')
legend(string(propA_vals))

%%
%Plot "benefit" of maintaining partner diversity 


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

% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N


%initialize results
leakiness_vals = [0 .5]; 
env_period_vals = [365/5 365/4 365/3 365/2 365 365*2 365*3 365*4 365*5 365*7];
propA_vals = 0:.1:1; 

results = nan(length(propA_vals), length(env_period_vals)); 

tspan = [1 8000]; 
count = 1;
count2 = 1;

figure(1)

d = 3
    difference_vals = [.6 .8 1]; 
    difference_val = difference_vals(4-d); 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;



for i = 1:length(env_period_vals)
    env_period = env_period_vals(i); 
    
    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol1 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, 1, 0, 0, 1, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res1 = deval(sol1, tspan(2)-env_period*3:tspan(2));

        sol2 = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, .5, .5, .5, .5, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res2 = deval(sol2, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res1(1,:)) - mean(final_res2(1,:));

            xticks([0:365:5000])
            xlabel('Days')
            
        
    end


end

figure(5)
pcolor(env_period_vals, propA_vals, results)
ylabel('PropA')
xlabel('Period')
h = colorbar
h.Label.String = 'Cost of Maintaining Diversity'






