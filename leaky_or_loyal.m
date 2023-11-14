%% Leaky or Loyal Model

% BLFS attempting to translate from Laura & Holly's model in R into MATLAB
% Nov 2023

% main System of Differential Equations

% Inputs:
% x, initial state space of system
        % x is size 1x7 
            % P, plant, photosynthetic biomass within tree
            % C, carbon belowground, Carbon pool from tree for fungal uptake
            % F1, fungus 1, symbiont biomass
            % F2, fungus 2, symbiont biomass
            % N, nitrogen aboveground, Nitrogen pool given to tree

% model parameters:
        % g, growth of plant proportional to Nitrogen pool
        % a, allocation of Carbon to Carbon pool
        % s, sesesence of tree biomass
        % l, loss of carbon pool to environment
        % r1_A, reward rate to fungus 1 in environment type A
        % r1_B, reward rate to fungus 1 in environment type B
        % r2_A, reward rate to fungus 2 in environment type A
        % r2_B, reward rate to fungus 2 in environment type B
        % e1, efficiency of fungus 1 carbon uptake
        % e2, efficiency of fungus 2 carbon uptake
        % m1, fungus 1 mortality
        % m2, fungus 2 mortality
        % d1, density dependent mortality of F1
        % d2, density dependent mortality of F2
        % u1_A, uptake of Nitrogen by fungus 1 in environment type A
        % u1_B, uptake of Nitrogen by fungus 1 in environment type B
        % u2_A, uptake of Nitrogen by fungus 2 in environment type A
        % u2_B, uptake of Nitrogen by fungus 1 in environment type B
        % mN, loss of Nitrogen from tree's stores
        % Ntot, total nitrogen = N + Ns 
            % Nitrogen aboveground + nitrogen in sediment

%envA = environment funciton of time t, eg envA = @(t) mod(t, 2);
 % envA(t) ==1 means environment is type A,
 % envA == 0 means environment is type B
        
% Output:  dxdt = derivatives of state variables at time t


function dxdt = leaky_or_loyal(t, x, g, a, s, l, r1_A, r1_B, r2_A, r2_B, e1, e2, m1, m2, d1, d2, u1_A, u1_B, u2_A, u2_B, mN, Ntot, envA)

%initialize output
dxdt = nan(size(x)); 

%sync rewards to environment
if envA(t) == 1
    r1 = r1_A;
    r2 = r2_A; 
    u1 = u1_A; 
    u2 = u2_A; 
elseif envA(t) == 0
    r1 = r1_B; 
    r2 = r2_B; 
    u1 = u1_B; 
    u2 = u2_B; 
end

dxdt(1) =  g.*x(5) - a.*x(1) - s.*x(1); % dP/dt = gN - aP - sP
dxdt(2) = a.*x(1) - l.*x(2) - (r1.*x(3) + r2.*x(4)).*x(2); % dC/dt = aP - lC - (r1F1 + r2F2)C
dxdt(3) = e1*r1.*x(3).*x(2) - m1.*x(3).*(1+d1*x(3)+d1*x(4)); % dF1/dt = e1*r1*F1*C - mF1(1+d1F1+d1F2)
dxdt(4) = e2*r2.*x(4).*x(2) - m2.*x(4).*(1+d2*x(4)+d2*x(3)); % dF2/dt = e2*r2*F2*C - mF2(1+d2F2+d2F1)
dxdt(5) = (u1.*x(3)+u2.*x(4)).*(Ntot-x(5)) - mN.*x(5); % dN/dt = (u1*F1 + u2*F2)(Ntot-N) - mN(N)


end
























