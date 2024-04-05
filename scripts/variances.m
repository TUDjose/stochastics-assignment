%% AE4304P - Stochastic Aerospace System Practical
% Analysis of simulated aircraft responses to atmospheric turbulence
% 
% José Bernardo Cunha (5216087)
%
% Modelling of symmetric aircraft dynamics and stability analysis

clc, clf, clear, close all;

model;
spectral_analysis;
close all;

format("shortG");


%% var.m function
y = lsim(A, B, C, D, u, t);
uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

varu = var(uhat);
varalpha = var(alpha);
vartheta = var(theta);
varqcV = var(qcV);
varnz = var(nz);

disp("matlab var.m")
disp([varu, varalpha, vartheta, varqcV, varnz])

y = lsim(Adamp, B, C, D, u, t);
uhat = y(:, 1);
alpha = y(:, 2);
theta = y(:, 3);
qcV = y(:, 4);

alphanz = alpha;
alphanz(length(alphanz) + 1) = 0;
nz = (V / g) * ((V / c) * qcV - diff(alphanz)/dt);
nz(length(nz)) = nz(length(nz) - 1);

varu = var(uhat);
varalpha = var(alpha);
vartheta = var(theta);
varqcV = var(qcV);
varnz = var(nz);

disp("matlab var.m damped")
disp([varu, varalpha, vartheta, varqcV, varnz])


%% analytic PSD
ana_var = zeros(1, 5);
for j=1:5
    for i=1:Nomega-1
        ana_var(j) = ana_var(j) + (w(i+1) - w(i)) * Sxx(i, j);
    end
end
disp("analytic var")
disp(ana_var/pi)

ana_var = zeros(1, 5);
for j=1:5
    for i=1:Nomega-1
        ana_var(j) = ana_var(j) + (w(i+1) - w(i)) * Sxx_damp(i, j);
    end
end
disp("analytic var damped")
disp(ana_var/pi)


%% periodograms
peri_var = zeros(1, 5);
for j=1:5
    for i=1:length(omega)-1
        peri_var(j) = peri_var(j) + (omega(i+1) - omega(i)) * periodograms(i, j);
    end
end
disp("exp var")
disp(peri_var/pi)

peri_var = zeros(1, 5);
for j=1:5
    for i=1:length(omega)-1
        peri_var(j) = peri_var(j) + (omega(i+1) - omega(i)) * periodograms_damp(i, j);
    end
end
disp("exp var damped")
disp(peri_var/pi)


%% smoothed periodograms
peri_var = zeros(1, 5);
for j=1:5
    F = smoothed(:, end) * 2 * pi;
    for i=1:length(F)-1
        peri_var(j) = peri_var(j) + (F(i+1) - F(i)) * smoothed(i, j);
    end
end
disp("exp smooth var")
disp(peri_var/pi)

peri_var = zeros(1, 5);
for j=1:5
    F = smoothed_damp(:, end) * 2 * pi;
    for i=1:length(F)-1
        peri_var(j) = peri_var(j) + (F(i+1) - F(i)) * smoothed_damp(i, j);
    end
end
disp("exp smooth var damped")
disp(peri_var/pi)


