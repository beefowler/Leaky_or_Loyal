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
e1 = 0.005; %efficiency of fungus 1 carbon uptake
e2 = 0.004; %efficiency of fungus 2 carbon uptake


rtot_vals = [.2 .8]; 
env_periods = [365 365*3 365*5 365*7];
colorval = cmocean('deep', 6); 
leakiness_vals = [-1:.2:1];

for rt = 1:2
    rtot = rtot_vals(rt) ;

    r1 = rtot / (e1/e2 + 1); 
    r2 = rtot - r1; 

    for k = 1:length(env_periods)
        env_period = env_periods(k);
        propA = .5;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;
        
        results = nan(1, length(leakiness_vals));
        fungus_ratio = results; 

        for i = 1:length(leakiness_vals)
            propleakiness = leakiness_vals(i);

            tspan = [1 8000];
            if tspan(2)-env_period*4 < 8000
                tspan = [1 env_period*8];
            end

            %simulate
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1*(1-propleakiness), r1*(1+propleakiness), r2*(1+propleakiness), r2*(1-propleakiness), e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:));

            results(i) = biomass_val ;
            fungus_ratio(i) = mean(final_res(3,:))./mean(final_res(4,:)); 

        
        subplot(2,2,rt)
        hold on
        plot(leakiness_vals, results, 'linewidth', 2, 'color', colorval(k+1, :));
        end


    end
end

    subplot(2,2,1)
    title('Low reward rate')
    ylabel('Tree biomass')
    xlabel('Leakiness')
    box on 

    legend({'1 year'; '3 years'; '5 years'; '7 years'}, 'Location', 'southwest')
    subplot(2,2,2)
    xlabel('Leakiness')
    title('High reward rate')
    box on 

%%

rtot_vals = [.2 .8]; 
env_periods = [365 365*3 365*5 365*7];
colorval = cmocean('deep', 6); 
e1 = 0.005; %efficiency of fungus 1 carbon uptake
e2 = 0.004; %efficiency of fungus 2 carbon uptake



for rt = 1:2
    rtot = rtot_vals(rt) ;

    r1 = rtot / (e1/e2 + 1); 
    r2 = rtot - r1; 

    leakiness_vals = -r2:.01:r2;


    for k = 1:length(env_periods)
        env_period = env_periods(k);
        propA = .5;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;
        
        results = nan(1, length(leakiness_vals));
        fungus_ratio = results; 

        for i = 1:length(leakiness_vals)
            leakiness = leakiness_vals(i);

            tspan = [1 8000];
            if tspan(2)-env_period*4 < 8000
                tspan = [1 env_period*8];
            end

            %simulate
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1-leakiness, r1+leakiness, r2+leakiness, r2-leakiness, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:));

            results(i) = biomass_val ;
            fungus_ratio(i) = mean(final_res(3,:))./mean(final_res(4,:)); 

        end

        subplot(2,2,rt)
        hold on
        plot(leakiness_vals, results, 'linewidth', 2, 'color', colorval(k+1, :));
    end
end

    subplot(2,2,1)
    title('Low reward rate')
    ylabel('Tree biomass')
    xlabel('Leakiness')
    box on 

    legend({'1 year'; '3 years'; '5 years'; '7 years'}, 'Location', 'southwest')
    subplot(2,2,2)
    xlabel('Leakiness')
    title('High reward rate')
    box on 

%% Now Mortality

figure 
rtot_vals = [.2 .8]; 
env_periods = [365 365*3 365*5 365*7];
colorval = cmocean('deep', 6); 
e1 = 0.005; %efficiency of fungus 1 carbon uptake
e2 = 0.005; %efficiency of fungus 2 carbon uptake

m1 = 0.005; 
m2 = 0.004; 

for rt = 1:2
    rtot = rtot_vals(rt) ;

    r1 = rtot / (m2/m1 + 1); 
    r2 = rtot - r1; 

    leakiness_vals = -r2:.01:r2;


    for k = 1:length(env_periods)
        env_period = env_periods(k);
        propA = .5;
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ;
        
        results = nan(1, length(leakiness_vals));
        fungus_ratio = results; 

        for i = 1:length(leakiness_vals)
            leakiness = leakiness_vals(i);

            tspan = [1 8000];
            if tspan(2)-env_period*4 < 8000
                tspan = [1 env_period*8];
            end

            %simulate
            sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1-leakiness, r1+leakiness, r2+leakiness, r2-leakiness, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
            final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
            biomass_val = mean(final_res(1,:));

            results(i) = biomass_val ;
            fungus_ratio(i) = mean(final_res(3,:))./mean(final_res(4,:)); 

        end

        subplot(2,2,rt)
        hold on
        plot(leakiness_vals, results, 'linewidth', 2, 'color', colorval(k+1, :));
    end
end

    subplot(2,2,1)
    title('Low reward rate')
    ylabel('Tree biomass')
    xlabel('Leakiness')
    box on 

    legend({'1 year'; '3 years'; '5 years'; '7 years'}, 'Location', 'southwest')
    subplot(2,2,2)
    xlabel('Leakiness')
    title('High reward rate')
    box on 



