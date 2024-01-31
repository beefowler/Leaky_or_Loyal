%Parameter space explore supplement

%What if fungi aren't symmetric in their other parameters 


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

tspan = [1 5000]; 

env_period = 365; 
propA = .5; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 



%% carbon conversion efficiency for each fungus


%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
m1 = 0.004; 
m2 = 0.004; 
mN = .1; 
Ntot = 10; 
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake

parameter1_vals = [0.001 0.005 0.01 0.05 0.1];
parameter2_vals = [0.001 0.005 0.01 0.05 0.1];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    e1 = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       e2 = parameter2_vals(j); 
       
        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 

       
        if  j == 2 & i == 1
            figure(2)
            subplot(1,2,1)
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
        end


    end

end

figure(3)

subplot(2,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
ylabel('m1')
set(gca, 'ydir', 'normal')


subplot(2,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
set(gca, 'ydir', 'normal')

colormap(cmocean('tempo'))

% carbon efficiency for each fungus -- bet-hedging strategy

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
m1 = 0.004; 
m2 = 0.004; 
mN = .1; 
Ntot = 10; 


parameter1_vals = [0.001 0.005 0.01 0.05 0.1];
parameter2_vals = [0.001 0.005 0.01 0.05 0.1];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    e1 = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       e2 = parameter2_vals(j); 
       
        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 

         if j == 2 & i == 1
            figure(2)
            subplot(1,2,2)
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
         end

    end

end

figure(3)
subplot(2,2,3)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('e2')
ylabel('e1')
set(gca, 'ydir', 'normal')


subplot(2,2,4)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('e2')
set(gca, 'ydir', 'normal')



colormap(cmocean('tempo'))


%% mortality for each fungus

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment

mN = .1; 
Ntot = 10; 
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake

parameter1_vals = [0.001 0.005 0.01 0.05 0.1];
parameter2_vals = [0.001 0.005 0.01 0.05 0.1];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    m1 = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       m2 = parameter2_vals(j); 
       
        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 

       
        if  j == 3 & i == 1
            figure(2)
            subplot(1,2,1)
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
        end


    end

end

figure(3)

subplot(2,2,1)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
ylabel('m1')
set(gca, 'ydir', 'normal')


subplot(2,2,2)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
set(gca, 'ydir', 'normal')

colormap(cmocean('tempo'))

% mortality for each fungus -- bet-hedging strategy

%defaults for all other parameters, 
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.01; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment

mN = .1; 
Ntot = 10; 
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake

parameter1_vals = [0.001 0.005 0.01 0.05 0.1];
parameter2_vals = [0.001 0.005 0.01 0.05 0.1];
results = nan(length(parameter1_vals), length(parameter2_vals)); 
var_results = results; 

for i = 1:length(parameter1_vals)
    m1 = parameter1_vals(i);

    for j = 1:length(parameter2_vals)
       m2 = parameter2_vals(j); 
       
        sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot/2, rtot/2, rtot/2, rtot/2, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(i,j) = mean(final_res(1,:)); 
        var_results(i,j) = std(final_res(1,:)); 

         if j == 3 & i == 1
            figure(2)
            subplot(1,2,2)
            b = plot(sol.x, sol.y([1 3 4], :));
            b(1).Color = [89,174,159]/255;
            b(2).Color =  [196/255 118/255 165/255];
            b(3).Color = [128, 180, 232]/255;

            b(1).LineWidth = 2;
            b(2).LineWidth = 2; 
            b(3).LineWidth = 2; 

            b(2).LineStyle = '-.';
            b(3).LineStyle = '--';

            legend({'Tree'; 'Fungus 1'; 'Fungus 2'})

            xticks([0:365:3000])
            xticklabels([0:3000/365])
            xlabel('Years')
            xlim([4*365 8*365])
         end

    end

end

figure(3)
subplot(2,2,3)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
ylabel('m1')
set(gca, 'ydir', 'normal')


subplot(2,2,4)
imagesc([1:length(parameter2_vals)], [1:length(parameter1_vals)], var_results)
xticks([1:length(parameter2_vals)])
xticklabels(parameter2_vals)
yticks([1:length(parameter1_vals)])
yticklabels(parameter1_vals)
xlabel('m2')
set(gca, 'ydir', 'normal')



colormap(cmocean('tempo'))
