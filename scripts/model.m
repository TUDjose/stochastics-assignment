%% AE4304P - Stochastic Aerospace System Practical
% Analysis of simulated aircraft responses to atmospheric turbulence
% 
% José Bernardo Cunha (5216087)
%
% Modelling of symmetric aircraft dynamics and stability analysis

clear;
clc;

%% aicraft constants
W = 44675;      % N
m = 4556;       % kg
S = 24.2;       % m^2
c = 2.022;      % m
b = 13.36;      % m
V = 51.4;       % m/s
h = 0;          % m -> LANDING CONFIGURATION
rho_air = 1.225;    % kg/m^2
lh = 5.5;      % m
muc = 76;
mub = 11;
KX2 = 0.012;
KZ2 = 0.037;
KZX = 0.002;
KY2 = 0.98;
xcg = 0.3 * c;  % m
g = 9.80665;


%% control derivatives
CX0 = 0;
CXu = -0.2173;
CXa = 0.4692;
CXq = 0;
CXde = 0;

CZ0 = -1.136;
CZu = -2.272;
CZa = -5.13;
CZadot = -1.405;
CZq = -3.84;
CZde = -0.6238;

Cmu = 0;
Cma = -0.4;
Cmadot = -3.615;
Cmq = -7.35;
Cmde  = -1.553;

CYb = -0.9896;
CYp = -0.087;
CYr = 0.43;
CYda = 0;
CYdr = 0.3037;

Clb = -0.0772;
Clp = -0.3415;
Clr = 0.283;
Clda = -0.2349;
Cldr = 0.0286;

Cnb = 0.1628;
Cnp = -0.0108;
Cnr = -0.1930;
Cnda = 0.0286;
Cndr = -0.1261;

Cmac = 0;
Cmh = 0;


%% turbulence parameters
Lg = 1500;      % m
sigma_wg = 2;   % m/s
sigma_ug = sigma_wg / V;   
sigma_ag = sigma_wg / V;   

CXug = CXu;
CXag = CXa;
CZug = CZu;
CZudotg = 2 * Cmac;
CZag = CZa;
CZadotg = CZadot  - CZq;
Cmug = Cmu;
Cmudotg = -2 * Cmh  * lh /c;
Cmag = Cma;
Cmadotg = Cmadot - Cmq;


%% full symmetric state space representation
xu = V * CXu  / (c * 2 * muc);
xa = V * CXa  / (c * 2 * muc);
xtheta = V * CZ0  / (c * 2 * muc);
xug = V * CXug  / (c * 2 * muc);
xag = V * CXag  / (c * 2 * muc);
xde = V * CXde  / (c * 2 * muc);

zu = V * CZu / (c * (2 * muc - CZadot));
za = V * CZa / (c * (2 * muc - CZadot));
ztheta =  V * -CX0 / (c * (2 * muc - CZadot));
zq = V * (2 * muc + CZq) / (c * (2 * muc - CZadot));
zug = V * CZug / (c * (2 * muc - CZadot));
zudotg = V * CZudotg / (c * (2 * muc - CZadot));
zag = V * CZag / (c * (2 * muc - CZadot));
zadotg = V * CZadotg / (c * (2 * muc - CZadot));
zde = V * CZde / (c * (2 * muc - CZadot));

mu = (V / c) * (Cmu + CZu * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
ma = (V / c) * (Cma + CZa * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
mtheta = (V / c) * (-CZ0 * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
mq = (V / c) * (Cmq + Cmadot * (2 * muc + CZq) / (2 * muc - CZadot)) / (2 * muc * KY2);
mug = (V / c) * (Cmug + CZug * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
mudotg = (V / c) * (Cmudotg + CZudotg * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
mag = (V / c) * (Cmag + CZag * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
madotg = (V / c) * (Cmadotg + CZadotg * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);
mde = (V / c) * (Cmde + CZde * Cmadot / (2 * muc - CZadot)) / (2 * muc * KY2);

% x =  [uhat, alpha, theta, qc/V, uhat_g, alpha_g, alphastar_g]
% u = [delta_e, w1, w3]

A = [xu xa xtheta 0 xug xag 0;
    zu za ztheta zq zug-zudotg*V*c/(Lg*V) zag zadotg*c/V;
    0 0 0 V/c 0 0 0;
    mu ma mtheta mq mug-mudotg*V*c/(Lg*V) mag madotg*c/V;
    0 0 0 0 -V/Lg 0 0;
    0 0 0 0 0 0 1;
    0 0 0 0 0 -(V/Lg)^2 -2*V/Lg];

B = [xde 0 0 ;
    zde zudotg*(c/V)*sigma_ug*sqrt(2*V/Lg) zadotg*(c/V)*sigma_ag*sqrt(3*V/Lg);
    0 0 0;
    mde mudotg*(c/V)*sigma_ug*sqrt(2*V/Lg) madotg*(c/V)*sigma_ag*sqrt(3*V/Lg);
    0 sigma_ug*sqrt(2*V/Lg) 0;
    0 0 sigma_ag*sqrt(3*V/Lg);
    0 0 (1-2*sqrt(3))*sigma_ag*sqrt((V/Lg)^3)];

C = eye(7,7);
D = zeros(7, 3);

Ktheta = -0.099;
K = [0 0 Ktheta 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0];
Adamp = A - B*K;

E = eig(Adamp);
lambda_sp = E(3);
zeta_sp = -real(lambda_sp) / abs(lambda_sp); 
% disp(zeta_sp)

sys = ss(Adamp, B, C, D);


%% add load factor to state space (w/ and wo/ pitch damper)
% nz = V/g * (thetadaot - alphadot)

nz_col =  (V / 9.80665) *  (A(:, 3) - A(:, 2));
Aload = [A nz_col;
    0 0 0 0 0 0 0 0];
Bload =[B; B(3, :) - B(2, :)];
Cload = eye(8,8);
Dload = zeros(8,3);

nz_col_damp = (V / 9.80665) *  (Adamp(:, 3) - Adamp(:, 2));
Aload_damp = [Adamp nz_col_damp;
                0 0 0 0 0 0 0 0];


%% plotting 
% disp(eig(Adamp))
% figure(3);
% pzmap(sys)
% 
% figure(1);
% step(sys);
% [y, t] = step(sys);
% y11 = y(1:2500, 4, 1);
% plot(t(1:2500), y11)
% xlabel('Time (s)')
% ylabel('Pitch rate (rad/s)')
% 
% figure(2);
% impulse(sys);
% [y, t] = impulse(sys);
% y14 = y(1:150, 4, 1);
% plot(t(1:150), y14)
% xlabel('Time (s)')
% ylabel('Pitch rate (rad/s)')





