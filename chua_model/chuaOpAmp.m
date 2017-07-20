function out = chuaOpAmp(istate,time,R,L)
    %%% istate  : vector of initial state [V1_init V2_init IL_init]
    %%% time    : total number of steps
    %%% R and L : resistance and optional inductance value
    %%% out     : timex3 matrix of all states for V1 V2 IL
    %%%
    %%% ex. run : out = chuaOpAmp([0.1 -0.1 0], 3000,1800);
    %%% plot    : plot3(out(:,1),out(:,2),-10 * out(:,3))
    
    % default inductance value
    if nargin < 4
        L = 1.8*10^(-2); % 18mH
    end
   
    % default step size
    step = 0.5*10^-5;
    h = step;
    
    % initializing output
    out = [istate;zeros(time,3)];
    state = istate;
    
    %%% 4th order Runga-Kutta %%%
    for i = 1:time
        k1 = h*Chua_ODE(state,R,L);
        k2 = h*Chua_ODE(state+k1/2,R,L);
        k3 = h*Chua_ODE(state+k2/2,R,L);
        k4 = h*Chua_ODE(state+k3,R,L);
        state = state + k1/6 + k2/3 + k3/3 + k4/6;
        out(i+1,:) = state;
    end
end


%%%%      Chua Equations        %%%%
function next_state = Chua_ODE(state,R,L)

    %%States to be estimated%%
    next_state = zeros(size(state));
    x = state(1);
    y = state(2);
    z = state(3);

    %%% circuit parameters %%%
    C1 = 10*10^(-9);   % 10  nF  
    C2 = 100*10^(-9);  % 100 nF
    r  = 10;           % 10 ohm
    G  = 1/R;
    
    % op amp rails
    vmax = 10.66;
    vmin = -10.16;

    %%% Chua diode parameters %%%
    %circuit values   theoretical values
    R1 = 3240;           %3300  ohms
    R2 = 21800;          %22000 ohms
    R3 = 21740;          %22000 ohms
    R4 = 2160;           %2200  ohms
    R5 = 219.7;            %220   ohms
    R6 = 221;            %221   ohms

    Pos_Esat = vmax * .57;
    Neg_Esat = -vmin * .72;

    % E1 and E2 define limit of piecewise functions
    Pos_E1 = R4/(R5+R4)*Pos_Esat;
    Pos_E2 = 1.2 * R1/(R2+R1)*Pos_Esat;

    Neg_E1 = R4/(R5+R4)*Neg_Esat;
    Neg_E2 = 1.1 * R1/(R2+R1)*Neg_Esat;

    % piecewise function slopes
    Ga = (-1/R1 - 1/R4);
    Gb = (1/R3 - 1/R4);
    Gc = .5 * (1/R6 + 1/R3);
    
    % piecewise intercepts
    Pos_b = 1.2 * Pos_E2*(Ga-Gb);
    Pos_c = Pos_E1*(Gc-Gb) + Pos_E2*(Gb-Ga) +.00007;
    Neg_b = Neg_E2*(Ga-Gb) + .00007;
    Neg_c = Neg_E1*(Gc-Gb) + Neg_E2*(Gb-Ga)-.00007;

    % multiplier for sigmoids
    mult = 1000;

    % current sigmoids
    I0 = (1./(1+exp(mult.*(x+Neg_E1)))) .* (x.*Gc + Neg_c);
    I1 = (1./(1+exp(mult.*(x+Neg_E2)))-1./(1+exp(mult.*(x+Neg_E1)))) .* (x.*Gb - Neg_b);
    I2 = (1./(1+exp(mult.*(x-Pos_E2)))-1./(1+exp(mult.*(x+Neg_E2)))) .* (x.*Ga-.00007);
    I3 = (1./(1+exp(mult.*(x-Pos_E1)))-1./(1+exp(mult.*(x-Pos_E2)))) .* (x.*Gb + Pos_b);
    I4 = (1./(1+exp(mult.*(-x+Pos_E1)))) .* (x.*Gc - Pos_c);

    % total current in Chua diode
    I = I0 + I1 + I2 + I3 + I4;
    
    % Chua's Circuit Equations
    next_state(1) = (1/C1)*(G*(y-x)-I);
    next_state(2) = (1/C2)*(G*(x-y)+z);
    next_state(3) = -(1/L)*z*r-(1/L)*y;
end
