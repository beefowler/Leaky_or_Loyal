

initial_conditions = [5]; 

for d = 1:length(initial_conditions)
  
[X2,FVAL2,~,~] = fmincon(@(x) negbiomass(x, initial_conditions(d)), .5, [],[],[],[],0,1); 

% 
% Xvec = [X1 X2 X3]; 
% Fvec = [FVAL1 FVAL2 FVAL3];
% [~,a] = min(Fvec);
% optimal_leakiness_values(d, 1) = Xvec(a); 
% optimal_leakiness_values(d, 2) = Fvec(a);

optimal_leakiness_values(d, 1) = X2; 
optimal_leakiness_values(d, 2) = FVAL2;

% go get simulation for this optimum 

g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d1_1 = 0.05;  %density dependent mortality effect of F1 on F1
d2_2 = 0.05;  %density dependent mortality effect of F2 on F2
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns
rtot = .2;

difference_val = 1; 
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;

x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = initial_conditions(d); %F1
x0(4) = 1; %F2
x0(5) = 1; %N

env_period = 365; 
tend = max(8000, env_period.*6);
tspan = [1 tend];
propA = 0.7; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

    leakiness = X2;
    r1_A = (1-leakiness).*rtot;
    r1_B = leakiness*rtot; 
    r2_A = leakiness*rtot; 
    r2_B = (1-leakiness)*rtot; 

   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   final_res = deval(sol, tspan(2)-env_period*3:tspan(2));
   %optimal_leakiness_sol = final_res; 
    optimal_leakiness_sol= deval(sol, tspan(1):tspan(end));

% Now Responsive solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, 0, 0, rtot, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   %responsive_sol = deval(sol, tspan(2)-env_period*3:tspan(2));
      responsive_sol = deval(sol, tspan(1):tspan(end));


% Now Loyal solution 
   sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, rtot, rtot, 0, 0, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
   %loyal_sol = deval(sol, tspan(2)-env_period*3:tspan(2));
   loyal_sol = deval(sol, tspan(1):tspan(end));


figure
subplot(1,3,1)
%b = plot(tspan(2)-env_period*3:tspan(2), loyal_sol);
b = plot(tspan(1):tspan(end), loyal_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title('Loyal')

subplot(1,3,2)
%b = plot(tspan(2)-env_period*3:tspan(2), responsive_sol);
b = plot(tspan(1):tspan(end), responsive_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title('Responsive')

subplot(1,3,3)
%b = plot(tspan(2)-env_period*3:tspan(2), optimal_leakiness_sol);
b = plot(tspan(1):tspan(end), optimal_leakiness_sol);
b(1).Color = [0 .9 .3];
b(2).Color = [0 .8 .7];
b(3).Color = 'r';
b(4).Color = 'b';
b(5).Color = [.5 0 .5];

b(2).LineStyle = '--';
b(5).LineStyle = '--';
title(['Optimally Leaky (' num2str(leakiness.*100) '%)'])

end

figure
subplot(1,2,1)
plot(initial_conditions, optimal_leakiness_values(:,1))
ylabel('Optimal leakiness')
subplot(1,2,2)
plot(initial_conditions, -optimal_leakiness_values(:,2))

ylabel('Optimal mean biomass')
xlabel('Initial Fungus 1 Population')




function B = negbiomass(leakiness, x0_F)
g = .1; %growth of plant proportional to Nitrogen pool
a = .04; %allocation of Carbon to Carbon pool
s = 0.00; %senesence of tree biomass
l = 0.005; %loss of carbon pool to environment
e1 = 0.05; %efficiency of fungus 1 carbon uptake
e2 = 0.05; %efficiency of fungus 2 carbon uptake
m1 = 0.001; %fungus 1 mortality
m2 = 0.001; %fungus 2 mortality
d2_1 = 0.01;  %density dependent mortality effect of F2 on F1
d1_2 = 0.01;  %density dependent mortality effect of F1 on F2
d1_1 = 0.05;  %density dependent mortality effect of F1 on F1
d2_2 = 0.05;  %density dependent mortality effect of F2 on F2
mN = .1; %loss of Nitrogen from tree's stores
Ntot = 10; %total nitrogen = N + Ns
rtot = .2;

difference_val = 1; 
u1_A = difference_val; %uptake of Nitrogen by fungus 1 in environment type A
u1_B = 1-difference_val; %uptake of Nitrogen by fungus 1 in environment type B
u2_A = 1-difference_val; %uptake of Nitrogen by fungus 2 in environment type A
u2_B = difference_val;

r1_A = (1-leakiness).*rtot;
r1_B = leakiness*rtot;
r2_A = leakiness*rtot;
r2_B = (1-leakiness)*rtot;

x0(1) = 1; %P
x0(2) = 1; %C
x0(3) = x0_F; %F1
x0(4) = 1; %F2
x0(5) = 1; %N

env_period = 365; 
tend = max(4000, env_period.*6);
tspan = [1 tend];
propA = 0.7; 
envA_treat = @(t) discretize(rem(t, env_period), [0 propA*env_period env_period]) == 1 ; 

sol = ode45(@(t, x) leaky_or_loyal_coexistence(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1_1, d2_1, d1_2, d2_2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA_treat), tspan, x0);
final_res = deval(sol, tspan(2)-env_period*3:tspan(2));

B = -mean(final_res(1,:)); %negative mean biomass


end

