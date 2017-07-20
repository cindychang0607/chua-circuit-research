function [stimpoint, jump] = Chuap(istate,timestep ,R ,dist,stimmag ,stimtime, dir)
% INPUT ARGUMENTS
% istate [vector]: initial conditions [xinit ; yinit ; zinit]
% timestep: step size in the model
% R : resistance value in model
% dist: distance on attractor from end of transient to perturbation
% stimag : pertubation magnitude
% stimtime: number of steps of perturbations
% dir [vector]: direction of pertubation in 1x3 vector
%
% OUTPUT ARGUMENTS
% stimpoint[vector]: coordinate of last perturbation
% jump: 1 if model tansitioned to other attractor, 0 otherwise 
%
%TEST
%[stimpoint jump] = Chuap([0;0.1;0], 1.2*10^-6, 1700, 4, .005 ,8,[0 0 1]);

state = istate;
substeps = 1;
dt = timestep;
h = dt/substeps;
distance = 0;
direction = [0 0 0];

laststate = state;
%%%% runs until distance is reached      
while distance < dist
    for j = 1:substeps
            k1 = h*Chua_ODE(state, R, direction);
            k2 = h*Chua_ODE(state+k1/2, R, direction);
            k3 = h*Chua_ODE(state+k2/2, R, direction);
            k4 = h*Chua_ODE(state+k3, R, direction);
            state = state+k1/6+k2/3+k3/3+k4/6;
            distance = distance + norm(state-laststate); % 3.3, 6.6, 9.9 ... multiples of 3.3
            laststate = state;
    end 
end

stimpoint = state;

%%%% perturbs model 
for kick = [0:stimmag*2/stimtime:stimmag stimmag:-2*stimmag/stimtime:0]
    direction = kick*dir;
    k1 = h*Chua_ODE(state, R, direction);
    k2 = h*Chua_ODE(state+k1/2, R, direction);
    k3 = h*Chua_ODE(state+k2/2, R, direction);
    k4 = h*Chua_ODE(state+k3, R, direction);
    state = state+k1/6+k2/3+k3/3+k4/6;
end

%%%% runs model after perturbation
xx = [state(1,:) zeros(1,2000)]; % stores V1 values 
for k = 1:2000
    k1 = h*Chua_ODE(state, R, direction);
    k2 = h*Chua_ODE(state+k1/2, R, direction);
    k3 = h*Chua_ODE(state+k2/2, R, direction);
    k4 = h*Chua_ODE(state+k3, R, direction);
    state = state+k1/6+k2/3+k3/3+k4/6;
    xx(k+1) = state(1,:);
end

jump = 0;
stdev = std(xx);
avg = mean(xx);
%%% mean of V1 values after perturbation
if -4 < avg && avg < -1 && 1.5 < stdev && stdev < 2.5
    jump = 1;
elseif .2 < avg && avg < 1 && 6 < stdev
    jump = 2;
end
end


%%%% Chua Model %%%%
function dx = Chua_ODE(state, R, kick)

% state varaibles with perturbation
dx = zeros(size(state));
x = state(1) + kick(1);
y = state(2) + kick(2);
z = state(3) + kick(3);

%%% Circuit Parameters %%%
C1 = 9.39*10^(-9);  %10    nF  
C2 = 96.8*10^(-9);  %100   nF 
r  = 10.3;          %10    ohms
R1 = 3239;          %3300  ohms
R2 = 18780;         %22000 ohms
R3 = 21740;         %22000 ohms
R4 = 2166;
R5 = 220.3; %220;
R6 = 221;   %220;%2200  ohms
L  = 1.431*10^(-2);  %33    mH %1.38mH
vmax = 9.829;
vmin = -8.656;

G = 1./R;
E1 = R4./(R5+R4)*vmax;
E2 = R1./(R2+R1)*vmax;

E11 = R4./(R5+R4)*(-vmin);
E22 = R1./(R2+R1)*(-vmin);
 
% % slopes
Ga = -1/R1 - 1/R4;
Gb =  1./R3 - 1/R4;
Gc =  1./R6 + 1./R3;

% shifts
b = E2.*(Ga-Gb);
c = E1.*(Gc-Gb) +  E2.*(Gb-Ga);

b1 = E22.*(Ga-Gb);
c1 = E11.*(Gc-Gb) +  E22.*(Gb-Ga);

%    Chua Equations
I0 = (1./(1+exp(1000*(x+E11)))) .* (x.*Gc + c1);
I1 = (1./(1+exp(1000*(x+E22)))-1./(1+exp(1000*(x+E11)))) .* (x.*Gb - b1);
I2 = (1./(1+exp(1000*(x-E2)))-1./(1+exp(1000*(x+E22)))) .* (x.*Ga);
I3 = (1./(1+exp(1000*(x-E1)))-1./(1+exp(1000*(x-E2)))) .* (x.*Gb + b);
I4 = (1./(1+exp(1000*(-x+E1)))) .* (x.*Gc - c);

%current
I = I0 + I1 + I2 + I3 + I4;
%   Chua's Circuit Equations
dx(1,:) = (1/C1).*(G.*(y-x)-I);
dx(2,:) = (1/C2).*(G.*(x-y)+z);
dx(3,:) = -(1./L).*z.*r-(1./L).*y;
end
