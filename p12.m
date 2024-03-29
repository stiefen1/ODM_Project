%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-483: Optimal Decision Making | Group Project, Exercise 1.2   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHORS | Bayane Benkhadda, Stephen Monnet, Bilel Hamrouni | 20.05.2022 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; yalmip('clear');clc;
%% Data
% The following data is define as a struct for each generator.
% You could define it in your own way if you want, 
% like vectors: GenCost=[15; 20; 15; 20; 30; 25];

% number of generators
NGen= 6; 
% time duration
T=24;
% cost
G1.cost = 15;
G2.cost = 20;
G3.cost = 15;
G4.cost = 20;
G5.cost = 30;
G6.cost = 25;
% capacity
G1.capacity = 10;
G2.capacity = 5;
G3.capacity = 10;
G4.capacity = 10;
G5.capacity = 20;
G6.capacity = 30;
% ramp-up
G1.rampup = 2;
G2.rampup = 5;
G3.rampup = 2;
G4.rampup = 5;
G5.rampup = 10;
G6.rampup = 5;
% ramp-down
G1.rampdown = 2;
G2.rampdown = 5;
G3.rampdown = 2;
G4.rampdown = 5;
G5.rampdown = 10;
G6.rampdown = 5;
% renewable energy source
r=[15.2;16.4;16.1;10.9;14.8;7.6;15.6;5.5;9.2;5.7;1.5;12.4;10.4;4.8;14.3;0.5;6.6;5.7;11.5;11.9;2.8;7.3;6.7;9.7];

% demand
d=[21.3;21.4;17.8;20.9;15.5;17.6;20.2;23.8;27.7;30.1;35.4;39.4;43.2;47.0;49.3;51.5;52.6;50.3;47.0;43.1;38.8;33.2;28.6;24.3];

% battery
battery.capacity = 20; % 20MWh
battery.ec = 0.95; % charging efficiency 
battery.ed = 0.92; % discharging efficiency

% Group to ease matricial computation
C = [G1.cost, G2.cost, G3.cost, G4.cost, G5.cost, G6.cost];
capacity = [G1.capacity, G2.capacity, G3.capacity, G4.capacity, G5.capacity, G6.capacity];
R_up = [G1.rampup, G2.rampup, G3.rampup, G4.rampup, G5.rampup, G6.rampup];
R_down = [G1.rampdown, G2.rampdown, G3.rampdown, G4.rampdown, G5.rampdown, G6.rampdown];

%% Economic Dispatch
% three elements
% Decision variables = Power deployed by each generator at each time-step
g = sdpvar(NGen, T);
% Charging decision
bc = sdpvar(1, T);
% Discharging decision
bd = sdpvar(1, T);
% Battery state
S = sdpvar(1, T);

con=[];%constraints initial
obj=0;%objective function initial

% objective function
obj = sum(C*g);

% constraints
for t = 1:1:T
    % Balance between production and consumption
    con = [con d(t) == sum(g(:,t)) + r(t) + bd(t) - bc(t) bd(t) >= 0 bc(t) >= 0];    
    
    for i = 1:1:NGen
        % Capacity constraints
        con = [con 0 <= g(i, t) g(i,t) <= capacity(i)];
        
        % Variation of power production
        if(t>1)
            con = [con g(i, t)-g(i, t-1) <= R_up(i) g(i, t-1)-g(i, t) <= R_down(i)];
        end
    end
    
    if(t>1)
        % Battery state
        con = [con S(t) == S(t-1) + battery.ec * bc(t-1) - bd(t-1)/battery.ed S(t) <= battery.capacity S(t) >= 0];
    end
end

% Enforce inital and final state of the battery to be 0
con = [con S(1) == 0 S(T) + battery.ec * bc(T) - bd(T)/battery.ed == 0];

%% define sdpsetting
ops=sdpsettings('solver','LINPROG');
sol=solvesdp(con,obj,ops);

% obtain the solutions and objective value
g_opt = value(g);
S_opt = value(S);
bc_opt = value(bc);
bd_opt = value(bd);
obj_opt = value(obj);

%% Verification of the solution
% Total produced power
power_production = sum(g_opt, 1)' + r;
time = [0:1:T-1];

% Check the production vs consumption equilibrium
figure;
plot(time, power_production, time, d, 'x', time, S_opt);
title("Power : production VS demand");
ylabel("Power [MW]");
xlabel("Time [h]");
legend("Production", "Demand", "Battery");
grid on;

g_der = zeros(NGen, T-1);
% Check the power variation constraint
for t = 2:1:T
    g_der(:, t-1) = g_opt(:, t) - g_opt(:, t-1);
end
 
% Select generator to check
Gen2check = 3;

figure;
plot(time(2:end), g_der(Gen2check, :), time, R_up(Gen2check)'.*[ones(1, T);-ones(1,T)], '--r');
ylim([-R_up(Gen2check)-1, R_up(Gen2check)+1]);
title(strcat("Power production of generator ", num2str(Gen2check)));
ylabel("Power [MW]");
xlabel("Time [h]");
grid on;
