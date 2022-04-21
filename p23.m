%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-483: Optimal Decision Making | Group Project, Exercise 2.2 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; yalmip('clear');clc;
%% Data
% number of generators
NGen= 6;
% time duration and time vector
T=24;
time = [0:1:T-1];

% generation cost
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
% start-up cost
G1.startup = 75;
G2.startup = 100;
G3.startup = 75;
G4.startup = 100;
G5.startup = 100;
G6.startup = 125;
% shutdown cost
G1.shutdown = 7.5;
G2.shutdown = 10;
G3.shutdown = 17.5;
G4.shutdown = 10;
G5.shutdown = 10;
G6.shutdown = 12.5;
% no-load cost
G1.noload = 10;
G2.noload = 5;
G3.noload = 10;
G4.noload = 10;
G5.noload = 10;
G6.noload = 10;
% initial state of generators x_0
G1.inital = 1;
G2.inital = 0;
G3.inital = 0;
G4.inital = 0;
G5.inital = 0;
G6.inital = 0;
% minimum up-time of generators
G1.minup = 3;
G2.minup = 3;
G3.minup = 3;
G4.minup = 3;
G5.minup = 3;
G6.minup = 3;
% minimum down-time of generators
G1.mindown = 2;
G2.mindown = 2;
G3.mindown = 2;
G4.mindown = 2;
G5.mindown = 2;
G6.mindown = 2;

% renewable energy source
r=[15.2;16.4;16.1;10.9;14.8;7.6;15.6;5.5;9.2;5.7;1.5;12.4;10.4;4.8;14.3;0.5;6.6;5.7;11.5;11.9;2.8;7.3;6.7;9.7];

% demand
d=[21.3;21.4;17.8;20.9;15.5;17.6;20.2;23.8;27.7;30.1;35.4;39.4;43.2;47.0;49.3;51.5;52.6;50.3;47.0;43.1;38.8;33.2;28.6;24.3];

% Define vectors to ease matricial computation
Cu = [G1.startup, G2.startup, G3.startup, G4.startup, G5.startup, G6.startup];
Cd = [G1.shutdown, G2.shutdown, G3.shutdown, G4.shutdown, G5.shutdown, G6.shutdown];
Cn = [G1.noload, G2.noload, G3.noload, G4.noload, G5.noload, G6.noload];
Cg = [G1.cost, G2.cost, G3.cost, G4.cost, G5.cost, G6.cost];
capacity = [G1.capacity, G2.capacity, G3.capacity, G4.capacity, G5.capacity, G6.capacity];
Tup = [G1.minup, G2.minup, G3.minup, G4.minup, G5.minup, G6.minup];
Tdn = [G1.mindown, G2.mindown, G3.mindown, G4.mindown, G5.mindown, G6.mindown];
x0 = [G1.inital, G2.inital, G3.inital, G4.inital, G5.inital, G6.inital]';

%% Unit commitment
x = binvar(NGen, T); % 1 = Gen up, 0 = Gen down
u = sdpvar(NGen, T); % 1 = Gen turned on, 0 = otherwise
v = sdpvar(NGen, T); % 1 = Gen turned off, 0 = otherwise
g = sdpvar(NGen, T);

% Initial state constraints
con = [x(:, 1) == x0 u(:, 1) == zeros(NGen, 1) v(:, 1) == zeros(NGen, 1)];

% objective function
obj = sum(Cu*u) + sum(Cd*v) + sum(Cn*x) + sum(Cg*g);

% Constraints over time
for t = 1:1:T
    % Balance between production and consumption
    con = [con d(t)==sum(g(:,t)) + r(t)];    
    
    for i = 1:1:NGen
        % Capacity constraints
        con = [con 0 <= g(i, t) g(i,t) <= capacity(i).*x(i,t)];
        
        % u and v are between 0 and 1
        con = [con u(i, t)>=0 v(i, t)>=0 u(i, t)<=1 v(i, t)<=1];
        
        if(t>1)
            % Constraint to be consistent in the state of the generator
            con = [con x(i,t-1)-x(i,t)+u(i,t) >= 0 x(i,t)-x(i,t-1)+v(i,t)>= 0];
            
            if(t<T)
                % Cons. 2f
                for tau = t+1:1:min(t+Tup(i), T)
                    con = [con x(i,t)-x(i,t-1) <= x(i,tau)];
                end
                
                % Cons. 2g
                for tau = t+1:1:min(t+Tdn(i), T)
                    con = [con x(i,t-1)-x(i,t) <= 1-x(i,tau)];
                end
            end
        end      
    end
end
%% define sdpsetting and solve
ops=sdpsettings('solver','GLPK', 'debug', 1);
sol=solvesdp(con,obj,ops);
solvertime=sol.solvertime;

% obtain the solutions and objective value
x_opt = value(x);
u_opt = value(u);
v_opt = value(v);
g_opt = value(g);

%% Plot results
% Plot generation of a choosed generator
Gen2plot = 2;

figure;
plot(time, g_opt(Gen2plot, :));
grid on;
title(strcat("Generator ", num2str(Gen2plot)));
xlabel("Time");
ylabel("Power Generation [MW]");

% Plot total production (generators + renewable) and demand
figure;
plot(time, sum(g_opt, 1), time, r', time, d, 'x');
grid on;
legend("Generator", "Renewable", "Demand");
title("Production VS Demand");
xlabel("Time");
ylabel("Power [MW]");
