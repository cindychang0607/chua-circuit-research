%%% Chua Diode./Op Amp Circuit IV Curve %%%

% op amp rails
vmax = 9;
vmin = -9;

%input voltage
step = (vmax-vmin)/10000; %step size
vin = vmin:step:vmax;     % Vin

%%% Chua diode parameters %%%
%circuit values   theoretical values
R1 = 3300;           %3300  ohms
R2 = 20000;          %22000 ohms
R3 = 20000;          %22000 ohms
R4 = 2200;           %2200  ohms
R5 = 220;            %220   ohms
R6 = 220;            %221   ohms

% Esaturation
Pos_Esat = vmax;
Neg_Esat = -vmin;

% E1 and E2 define limit of piecewise functions
Pos_E1 = R4/(R5+R4)*Pos_Esat;
Pos_E2 = R1/(R2+R1)*Pos_Esat;

Neg_E1 = R4/(R5+R4)*Neg_Esat;
Neg_E2 = R1/(R2+R1)*Neg_Esat;


% piecewise function slopes
Ga = (-1/R1 - 1/R4);
Gb = (1/R3 - 1/R4);
Gc = (1/R6 + 1/R3);

% piecewise intercepts
Pos_b = Pos_E2*(Ga-Gb);
Pos_c = Pos_E1*(Gc-Gb) + Pos_E2*(Gb-Ga);
Neg_b = Neg_E2*(Ga-Gb);
Neg_c = Neg_E1*(Gc-Gb) + Neg_E2*(Gb-Ga);

% multiplier for sigmoids
mult = 1000;

% current sigmoids
I0 = (1./(1+exp(mult.*(vin+Neg_E1)))) .* (vin.*Gc + Neg_c);
I1 = (1./(1+exp(mult.*(vin+Neg_E2)))-1./(1+exp(mult.*(vin+Neg_E1)))) .* (vin.*Gb - Neg_b);
I2 = (1./(1+exp(mult.*(vin-Pos_E2)))-1./(1+exp(mult.*(vin+Neg_E2)))) .* (vin.*Ga);
I3 = (1./(1+exp(mult.*(vin-Pos_E1)))-1./(1+exp(mult.*(vin-Pos_E2)))) .* (vin.*Gb + Pos_b);
I4 = (1./(1+exp(mult.*(-vin+Pos_E1)))) .* (vin.*Gc - Pos_c);

% total current in Chua diode
I = (I0 + I1 + I2 + I3 + I4);

plot(vin, I*10);
grid on;

