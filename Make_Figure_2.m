%% Leakiness in different environments 
% now vary reward strategy in fixed environmental regime
clear all
figure

% set parameter values
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.01; %efficiency of fungus 1 carbon uptake
e2 = 0.01; %efficiency of fungus 2 carbon uptake
m1 = 0.01; %fungus 1 mortality
m2 = 0.01; %fungus 2 mortality
d1 = 0.1; %density dependent mortality of F1
d2 = 0.1; %density dependent mortality of F2
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns


% Initial conditions
x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = 1; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

%initialize results
leakiness_vals = 0:.05:1; 
env_period_vals = [180 365 365*2 365*3];
propA = .5; 
results = nan(length(env_period_vals), length(leakiness_vals)); 
extinctionresults = nan(length(env_period_vals), length(leakiness_vals)); 

tspan = [1 8000]; 

    difference_val = 1; 
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;


    coloropts = viridis(10); 

    count = -1; 
for rtot = [.2 .8]; 
       count = count+2

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = rtot.*(1-leakiness);
    r1_B = rtot.*leakiness; 
    r2_A = rtot.*leakiness; 
    r2_B = rtot.*(1-leakiness); 

    for j = 1:length(env_period_vals)

        env_period = env_period_vals(j);
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));
    
        %check if one of the funguses is going extinct 
        tpoint_1 = env_period*2;
        tpoint_2 = env_period*3;
        if mean(final_res(3,:))>mean(final_res(4,:))
            if final_res(4,tpoint_2)<final_res(4,tpoint_1)
                extinctionresults(j,i) = 1; %1 is extinct, 0 is not extinct
            else
                extinctionresults(j,i) = 0;
            end
        else 
             if final_res(3,tpoint_2)<final_res(3,tpoint_1)
                extinctionresults(j,i) = 1; %1 is extinct, 0 is not extinct
            else
                extinctionresults(j,i) = 0;
             end
        end

    end


end

subplot(2,2,count)
plot(leakiness_vals, results(1,:), 'linewidth', 2, 'color', coloropts(10,:))
hold on 
plot(leakiness_vals, results(2,:),'linewidth', 2, 'color',  coloropts(8,:))
plot(leakiness_vals, results(3,:),'linewidth', 2, 'color',  coloropts(6,:))
plot(leakiness_vals, results(4,:),'linewidth', 2, 'color',  coloropts(4,:))

xlabel('Leakiness')
%ylim([0 16])

%scatter(leakiness_vals, results(1,:), 20, extinctionresults(1,:), 'filled', 'MarkerEdgeColor', 'k')
%scatter(leakiness_vals, results(2,:), 20, extinctionresults(2,:), 'filled', 'MarkerEdgeColor', 'k')
%scatter(leakiness_vals, results(3,:), 20, extinctionresults(3,:), 'filled', 'MarkerEdgeColor', 'k')
%scatter(leakiness_vals, results(4,:), 20, extinctionresults(4,:), 'filled', 'MarkerEdgeColor', 'k')
colormap('bone')
ylabel('Average tree biomass')

end

subplot(2,2,1)
ylabel('Average tree biomass')
title('Even environment, variable period')
legend({'6 months'; '1 year'; '2 years'; '3 years'}); 

subplot(2,2,2)
hold on 
title('Annual period, variable evenness')



%initialize results
env_period = 365;
propA_vals = .5:.1:1; 

results = nan(length(propA_vals), length(leakiness_vals)); 

count = 1;
count2 = 1;

rtot = 0.2;

    difference_val = 1
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = rtot.*(1-leakiness);
    r1_B = rtot.*leakiness; 
    r2_A = rtot.*leakiness; 
    r2_B = rtot.*(1-leakiness); 

    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

         % if propA == 0.9
         %     if leakiness == 0
         %     subplot(2,6,4)
         %     elseif leakiness == leakiness_vals(7)
         %     subplot(2,6,5)
         %     elseif leakiness == 0.5
         %     subplot(2,6,6)
         %     elseif leakiness == 0.55
         %     subplot(2,6,10)
         %     elseif leakiness == 0.8
         %     subplot(2,6,11)
         %     elseif leakiness == 1 
         %     subplot(2,6,12)
         %     else
         %         continue
         %     end
         %    b = plot(sol.x, sol.y);
         %    b(1).Color = [0 .9 .3];
         %    b(2).Color = [0 .8 .7];
         %    b(3).Color = 'r';
         %    b(4).Color = 'b';
         %    b(5).Color = [.5 0 .5];
         % 
         %    b(2).LineStyle = '--';
         %    b(5).LineStyle = '--';
         % 
         %    legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
         %    title(['leakiness = ' num2str(leakiness)])
         % 
         %    xticks([0:365:8000])
         %    xlabel('Days')
         % 
         % end
    end


end

subplot(2,2,2)
plot(leakiness_vals, results(1,:), 'linewidth', 2, 'color', coloropts(10,:))
hold on 
plot(leakiness_vals, results(2,:),'linewidth', 2, 'color',  coloropts(8,:))
plot(leakiness_vals, results(3,:),'linewidth', 2, 'color',  coloropts(7,:))
plot(leakiness_vals, results(4,:),'linewidth', 2, 'color',  coloropts(6,:))
plot(leakiness_vals, results(5,:),'linewidth', 2, 'color',  coloropts(4,:))
plot(leakiness_vals, results(6,:),'linewidth', 2, 'color',  coloropts(3,:))

%ylabel('Average tree biomass')
xlabel('Leakiness')
legend(string(propA_vals))

set(findall(gcf,'-property','FontSize'),'FontSize',13)

subplot(2,2,2)
ylim([0 25])
subplot(2,2,1)
ylim([0 25])
subplot(2,2,3)
ylim([0 25])


%%
count = 1;
count2 = 1;

rtot = 0.8;

    difference_val = 1
    u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
    u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
    u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
    u2_B = difference_val;

for i = 1:length(leakiness_vals)
    leakiness = leakiness_vals(i);

    r1_A = rtot.*(1-leakiness);
    r1_B = rtot.*leakiness; 
    r2_A = rtot.*leakiness; 
    r2_B = rtot.*(1-leakiness); 

    for j = 1:length(propA_vals)
        propA = propA_vals(j); 
        envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

        sol = ode45(@(t, x) leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
        final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
        results(j,i) = mean(final_res(1,:));

         % if propA == 0.9
         %     if leakiness == 0
         %     subplot(2,6,4)
         %     elseif leakiness == leakiness_vals(7)
         %     subplot(2,6,5)
         %     elseif leakiness == 0.5
         %     subplot(2,6,6)
         %     elseif leakiness == 0.55
         %     subplot(2,6,10)
         %     elseif leakiness == 0.8
         %     subplot(2,6,11)
         %     elseif leakiness == 1 
         %     subplot(2,6,12)
         %     else
         %         continue
         %     end
         %    b = plot(sol.x, sol.y);
         %    b(1).Color = [0 .9 .3];
         %    b(2).Color = [0 .8 .7];
         %    b(3).Color = 'r';
         %    b(4).Color = 'b';
         %    b(5).Color = [.5 0 .5];
         % 
         %    b(2).LineStyle = '--';
         %    b(5).LineStyle = '--';
         % 
         %    legend({'Tree'; 'Carbon allocation'; 'Fungus 1'; 'Fungus 2'; 'Nitrogen in Tree'})
         %    title(['leakiness = ' num2str(leakiness)])
         % 
         %    xticks([0:365:8000])
         %    xlabel('Days')
         % 
         % end
    end


end

subplot(2,2,4)
plot(leakiness_vals, results(1,:), 'linewidth', 2, 'color', coloropts(10,:))
hold on 
plot(leakiness_vals, results(2,:),'linewidth', 2, 'color',  coloropts(8,:))
plot(leakiness_vals, results(3,:),'linewidth', 2, 'color',  coloropts(7,:))
plot(leakiness_vals, results(4,:),'linewidth', 2, 'color',  coloropts(6,:))
plot(leakiness_vals, results(5,:),'linewidth', 2, 'color',  coloropts(4,:))
plot(leakiness_vals, results(6,:),'linewidth', 2, 'color',  coloropts(3,:))

ylabel('Average tree biomass')
xlabel('Leakiness')
legend(string(propA_vals))

set(findall(gcf,'-property','FontSize'),'FontSize',13)


subplot(2,2,4)
ylim([0 25])

