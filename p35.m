%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPFL | MGT-483: Optimal Decision Making | Group Project, Exercise 3.5 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; %close all; 
yalmip('clear');clc;
%% Data
% The following data is define as a struct for each generator.
% You could define it in your own way if you want, 
% like vectors: GenCost=[15; 20; 15; 20; 30; 25];

% number of generators
NGen= 6;
% time duration
T=24;
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
r_bar=[15.2;16.4;16.1;10.9;14.8;7.6;15.6;5.5;9.2;5.7;1.5;12.4;10.4;4.8;14.3;0.5;6.6;5.7;11.5;11.9;2.8;7.3;6.7;9.7];
r_hat = 0.6*ones(T,1);
%r_hat = 0.3*ones(T,1);
capacity = [G1.capacity, G2.capacity, G3.capacity, G4.capacity, G5.capacity, G6.capacity];

% demand
d=[21.3;21.4;17.8;20.9;15.5;17.6;20.2;23.8;27.7;30.1;35.4;39.4;43.2;47.0;49.3;51.5;52.6;50.3;47.0;43.1;38.8;33.2;28.6;24.3];

% Generation cost
Cg = [G1.cost, G2.cost, G3.cost, G4.cost, G5.cost, G6.cost];
Cu = [G1.startup, G2.startup, G3.startup, G4.startup, G5.startup, G6.startup];
Cd = [G1.shutdown, G2.shutdown, G3.shutdown, G4.shutdown, G5.shutdown, G6.shutdown];
Cn = [G1.noload, G2.noload, G3.noload, G4.noload, G5.noload, G6.noload];
Tup = [G1.minup, G2.minup, G3.minup, G4.minup, G5.minup, G6.minup];
Tdn = [G1.mindown, G2.mindown, G3.mindown, G4.mindown, G5.mindown, G6.mindown];

%% Robust Unit Commitment
% Hints: define binary variables (vector) in yalmip:  
% a = binvar(N,M) with dimension N*M 
x = binvar(T, NGen, 'full'); % 1 = Gen up, 0 = Gen down
u = binvar(T, NGen, 'full'); % 1 = Gen turned on, 0 = otherwise
v = binvar(T, NGen, 'full'); % 1 = Gen turned off, 0 = otherwise

% Pour chaque timestep pour chaque générateur
a = sdpvar(T, NGen, 'full');
b = sdpvar(T-1, NGen, 'full');
%b = zeros(T-1, NGen, 'like', sdpvar);
c = sdpvar(T, NGen, 'full');

% Dualized problem
y = sdpvar(2*T, 1);
z = zeros(T, 1, 'like', sdpvar);
add_cost = 0; % Cost that is excluded from the r \in R

for t=1:1:T-1
    z(t) = Cg * (a(t, :) + b(t, :))';
    add_cost = add_cost + Cg * c(t, :)';
end
z(T) = Cg * a(T, :)';

w = [r_bar + r_hat; r_bar + r_hat];
A = [eye(T); -eye(T)]';

% Constraint 3c
P = sdpvar(4, T-1, NGen);
Q = sdpvar(4, T-1, NGen);

P0 = sdpvar(2, NGen);
Q0 = sdpvar(2, NGen);

R = [r_bar(1:end-1) + r_hat(1:end-1), r_bar(1:end-1) + r_hat(1:end-1),...
    r_bar(2:end) + r_hat(2:end), r_bar(2:end) + r_hat(2:end)];
R0 = [r_bar(1) + r_hat(1), r_bar(1) + r_hat(1)];

C = [1 -1 0 0;...
    0 0 1 -1];
B = sdpvar(2, T-1, NGen);

for i=1:NGen
    B(1, :, i) = b(:, i)';
    B(2, :, i) = a(2:end, i)';
end


con=[];%constraints initial
obj=0;%objective function initial

%% CONSTRAINTS
% Equilibrium between production and demand
for t=1:1:T
    % Equilibrium   
    %con = [con abs(a(t, :)) <= x(t, :)*1e3];
    con = [con sum(a(t, :)) == -1];
    con = [con sum(c(t, :)) == d(t)];
    %con = [con abs(c(t, :)) <= x(t, :)*1e3];
    
    if t < T
        con = [con sum(b(t, :)) == 0];
        %con = [con abs(b(t, :)) <= x(t, :)*1e3];
    end
end

% Constraints from previous ex.23

for t = 1:1:T    
    for i = 1:1:NGen
        if(t>1)
            % Constraint to be consistent in the state of the generator
            con = [con x(t-1, i)-x(t, i)+u(t, i) >= 0 x(t, i)-x(t-1, i)+v(t, i)>= 0];
            
            if(t<T)
                % Cons. 2f
                for tau = t+1:1:min(t+Tup(i), T)
                    con = [con x(t, i)-x(t-1, i) <= x(tau, i)];
                end
                
                % Cons. 2g
                for tau = t+1:1:min(t+Tdn(i), T)
                    con = [con x(t-1, i)-x(t, i) <= 1-x(tau, i)];
                end
            end
        end      
    end
end

% r is in an uncertainty box
con = [con A*y == z];
con = [con y >= 0];

% Constraint 3c
for i=1:NGen
    % Constraints for t = 0
    con = [con [1 -1]*P0(:, i) == a(1, i)];
    con = [con R0*P0(:, i) >= -c(1, i)];
    
    con = [con [1 -1]*Q0(:, i) == -a(1, i)];
    con = [con R0*Q0(:, i) >= c(1, i) - capacity(i)*x(1,i)];
    
    % Constraints for t > 0
    con = [con C*P(:, :, i)==B(:, :, i)];
    diagRP = diag(R*P(:, :, i));
    con = [con diagRP >= -c(2:end, i)];
    
    con = [con C*Q(:, :, i)==-B(:, :, i)];
    diagRQ = diag(R*Q(:, :, i));
    con = [con diagRQ >= c(2:end, i) - capacity(i)*x(2:end, i)];

end
con = [con P <= 0];
con = [con Q <= 0];
con = [con P0 <= 0];
con = [con Q0 <= 0];

% Enforce first generator to be ON at t = 1
con = [con x(1, 1) == 1];

%% Objective function
for t=1:1:T
    for i=1:1:NGen
        obj = obj + Cu(i) * u(t, i) + Cd(i) * v(t, i) + Cn(i) * x(t, i);
    end
end

obj = obj + w'*y + add_cost;

%% define sdpsetting
ops=sdpsettings('solver','MOSEK');%, 'debug', 1);
sol=solvesdp(con,obj,ops);

% obtain the solutions and objective value
a_opt = value(a);
b_opt = value(b);
c_opt = value(c);

%% Plot results
prod = zeros(T, NGen);
r = r_bar + 2*(rand(T,1)-0.5).*r_hat; % r varies between r_bar-r_hat and r_bar+r_hat
for i = 1:1:NGen
    for t=2:1:T
        prod(t, i) = x(t, i)*(a_opt(t,i)*r(t) + b_opt(t-1,i)*r(t-1) + c_opt(t, i));
    end
    prod(1, i) =  x(1, i)*(a_opt(1,i)*r(1) + c_opt(1, i));
end
time = [1:1:T];
prod_tot = sum(prod, 2); % total production at each time-step

figure;
plot(time, prod_tot + r, '--r');
hold on;
plot(time, d, 'x', 'Color', 'b');
title("Equilibrium between production and demand");
ylabel("Power [MW]");
xlabel("Time [h]");
legend("Production", "Demand");
grid on;


%% Plot generation of a choosed generator
close gcf
Gen2plot = 1;

figure;
plot(time, prod(:, Gen2plot));
hold on;
yline([0 capacity(Gen2plot)], '--r');
grid on;
title(strcat("Generator ", num2str(Gen2plot)));
xlabel("Time [h]");
ylabel("Power Generation [MW]");
ylim([-0.5 capacity(Gen2plot)+0.5]);

%% Plot all the generators
figure;

for i=4:1:6
    subplot(3, 1, i-3);
    plot(time, prod(:, i));
    hold on;
    yline([0 capacity(i)], '--r');
    grid on;
    title(strcat("Generator ", num2str(i)));
    xlabel("Time [h]");
    ylabel("Power [MW]");
    ylim([-0.5 capacity(i)+0.5]);
end
