%% AE4304P - Stochastic Aerospace System Practical
% Analysis of simulated aircraft responses to atmospheric turbulence
% 
% José Bernardo Cunha (5216087)
%
% Modelling of symmetric aircraft dynamics and stability analysis

clc, clf, clear, close all;

model;
close all;

%% analytic PSD
Nomega = 300;
w = logspace(-2, 2, Nomega);

Cu =     [1 0 0 0 0 0 0];
Calpha = [0 1 0 0 0 0 0];
Ctheta = [0 0 1 0 0 0 0];
Cq =     [0 0 0 1 0 0 0];
Dr =      [0 0 0];

Calphadot = A(2, :);
Dalphadot = B(2, :);

mag = bode(A, B, Cu, Dr, 3, w); Suu = mag.*mag;
mag = bode(A, B, Calpha, Dr, 3, w); Saa = mag.*mag;
mag = bode(A, B, Ctheta, Dr, 3, w); Stt = mag.*mag;
mag = bode(A, B, Cq, Dr, 3, w); Sqq = mag.*mag;

sys = ss(A, B(:, 3), Calphadot, Dalphadot(:, 3));
Hadotw = freqresp(sys, w);
sys = ss(A, B(:, 3), Cq, Dr(:, 3));
Hqw = freqresp(sys, w);

Hnzw = (V / g) * ((V / c)*Hqw - Hadotw);
mag = abs(squeeze(Hnzw)); Snznz = mag.*mag;


Sxx = [Suu Saa Stt Sqq Snznz];
labels = ["$S_{uu} \frac{rad^2}{Hz}$"
    "$S_{\alpha\alpha} \frac{rad^2}{Hz}$"
    "$S_{\theta\theta} \frac{rad^2}{Hz}$"
    "$S_{qq} \frac{rad^2}{Hz}$"
    "$S_{nznz} \frac{rad^2}{Hz}$"];

labels_hat = ["$\hat{S}_{uu} \frac{rad^2}{Hz}$"
    "$\hat{S}_{\alpha\alpha} \frac{rad^2}{Hz}$"
    "$\hat{S}_{\theta\theta} \frac{rad^2}{Hz}$"
    "$\hat{S}_{qq} \frac{rad^2}{Hz}$"
    "$\hat{S}_{nznz} \frac{rad^2}{Hz}$"];

figure(1)
for i = 1:5
    subplot(3, 2, i);
    loglog(w, Sxx(:, i));
    xlabel('log $\omega$','Interpreter','latex');
    ylabel(labels(i),'Interpreter','latex')
end

mag = bode(Adamp, B, Cu, Dr, 3, w); Suu = mag.*mag;
mag = bode(Adamp, B, Calpha, Dr, 3, w); Saa = mag.*mag;
mag = bode(Adamp, B, Ctheta, Dr, 3, w); Stt = mag.*mag;
mag = bode(Adamp, B, Cq, Dr, 3, w); Sqq = mag.*mag;

Calphadot = Adamp(2, :);
Dalphadot = B(2, :);

sys = ss(Adamp, B(:, 3), Calphadot, Dalphadot(:, 3));
Hadotw = freqresp(sys, w);
sys = ss(Adamp, B(:, 3), Cq, Dr(:, 3));
Hqw = freqresp(sys, w);

Hnzw = (V / g) * ((V / c)*Hqw - Hadotw);
mag = abs(squeeze(Hnzw)); Snznz = mag.*mag;

Sxx_damp = [Suu Saa Stt Sqq Snznz];

figure(2)
for i = 1:5
    subplot(3, 2, i);
    loglog(w, Sxx_damp(:, i), w, Sxx(:, i), '--');
    xlabel('log $\omega$','Interpreter','latex');
    ylabel(labels(i),'Interpreter','latex')
end


%% experimental PSD method: fft
dt = 0.01;
fs = 1 / dt;
T = 60;
t = [0:dt:T];
N = length(t);

nn = zeros(1, N);
w3 = randn(1, N) / sqrt(dt);

u = [nn' nn' w3'];      % vertical turbulence and no elevator deflection

y = lsim(A, B, C, D, u, t);
uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

U = dt * fft(uhat);
ALPHA = dt * fft(alpha);
THETA = dt * fft(theta);
QCV = dt * fft(qcV);
NZ = dt * fft(nz);

Pu = (1/T) * (U .* conj(U));
Palpha = (1/T) * (ALPHA .* conj(ALPHA));
Ptheta = (1/T) * (THETA .* conj(THETA));
PqcV = (1/T) * (QCV .* conj(QCV));
Pnz = (1/T) * (NZ .* conj(NZ));

periodograms = [Pu Palpha Ptheta PqcV Pnz];

omega = 2*pi*fs*(0:(N/2)-1)/N;

figure(3);
for i = 1:5
    subplot(3,2,i);
    loglog(w, Sxx(:, i), '--', omega, periodograms(1 : round (N/2) - 1, i));
    ylim([10^-10, 1]);
    xlabel("log $\omega$", 'Interpreter', 'latex');
    ylabel(labels_hat(i), 'Interpreter', 'latex');
end


y = lsim(Adamp, B, C, D, u, t);
uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

U = dt * fft(uhat);
ALPHA = dt * fft(alpha);
THETA = dt * fft(theta);
QCV = dt * fft(qcV);
NZ = dt * fft(nz);

Pu = (1/T) * (U .* conj(U));
Palpha = (1/T) * (ALPHA .* conj(ALPHA));
Ptheta = (1/T) * (THETA .* conj(THETA));
PqcV = (1/T) * (QCV .* conj(QCV));
Pnz = (1/T) * (NZ .* conj(NZ));

periodograms_damp = [Pu Palpha Ptheta PqcV Pnz];

omega = 2*pi*fs*(0:(N/2)-1)/N;

figure(4);
for i = 1:5
    subplot(3,2,i);
    loglog(w, Sxx_damp(:, i), '--', omega, periodograms_damp(1 : round (N/2) - 1, i));
    ylim([10^-10, 1]);
    xlabel("log $\omega$", 'Interpreter', 'latex');
    ylabel(labels_hat(i), 'Interpreter', 'latex');
end


%% experimental PSD method: pwelch
y = lsim(A, B, C, D, u, t); 

uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

[puu, f] = pwelch(uhat, [], [], [], fs);
[paa, f] = pwelch(alpha, [], [], [], fs);
[ptt, f] = pwelch(theta, [], [], [], fs);
[pqq, f] = pwelch(qcV, [], [], [], fs);
[pnn, f] = pwelch(nz, [], [], [], fs);

smoothed = [puu, paa, ptt, pqq, pnn, f];

figure(5);
for i=1:5
    subplot(3,2,i);
    loglog(w, Sxx(:, i), '--', 2*pi*f,  smoothed(:, i));
    ylim([10^-10, 1]);
    xlabel("log $\omega$", 'Interpreter', 'latex');
    ylabel(labels_hat(i), 'Interpreter', 'latex');
end

y = lsim(Adamp, B, C, D, u, t); 

uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

[puu, f] = pwelch(uhat, [], [], [], fs);
[paa, f] = pwelch(alpha, [], [], [], fs);
[ptt, f] = pwelch(theta, [], [], [], fs);
[pqq, f] = pwelch(qcV, [], [], [], fs);
[pnn, f] = pwelch(nz, [], [], [], fs);

smoothed_damp = [puu, paa, ptt, pqq, pnn, f];

figure(6);
for i=1:5
    subplot(3,2,i);
    loglog(w, Sxx_damp(:, i), '--', 2*pi *f,  smoothed_damp(:, i));
    ylim([10^-10, 1]);
    xlabel("log $\omega$", 'Interpreter', 'latex');
    ylabel(labels_hat(i), 'Interpreter', 'latex');
end


