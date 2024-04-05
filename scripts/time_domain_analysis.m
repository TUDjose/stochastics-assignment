%% AE4304P - Stochastic Aerospace System Practical
% Analysis of simulated aircraft responses to atmospheric turbulence
% 
% José Bernardo Cunha (5216087)
%
% Modelling of symmetric aircraft dynamics and stability analysis

clc, clf, clear, close all;

model;
close all;

%% generate turbulence input vector
dt = 0.005;
T = 60;
t = [0:dt:T];
N = length(t);

nn = zeros(1, N);
w3 = sigma_wg .* randn(1, N) / sqrt(dt);

u = [nn' nn' w3'];      % vertical turbulence and no elevator deflection


%% state response
y = lsim(A, B, C, D, u, t);
y2 = lsim(Adamp, B, C, D, u, t);

%% plotting

figure(1);
subplot ( 5 , 1 , 1 ) ;
plot ( t , y ( : , 1 ));
xlabel('Time (s)'); ylabel('$u/V$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 2 ) ;
plot ( t , y ( : , 2 ) * (180/pi))  
xlabel('Time (s)'); ylabel('$\alpha$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 3 ) ;
plot ( t , y ( : , 3 )* (180/pi))
xlabel('Time (s)'); ylabel('$\theta$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 4 ) ;
plot ( t , y ( : , 4 )* (180/pi))
xlabel('Time (s)'); ylabel('$qc/V$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 5 ) ;
nz = (V / g) * ((V/c)*y(:, 4)' - A(2, :) * y');
plot(t, nz);
xlabel('Time (s)'); ylabel('$n_z$ (-)','Interpreter','latex')


figure(2);
subplot ( 5 , 1 , 1 ) ;
plot ( t , y2 ( : , 1 )', t , y ( : , 1 ), '--');
xlabel('Time (s)'); ylabel('$u/V$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 2 ) ;
plot ( t , y2 ( : , 2 ) * (180/pi),  t , y ( : , 2 ) * (180/pi), '--')  
xlabel('Time (s)'); ylabel('$\alpha$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 3 ) ;
plot ( t , y2 ( : , 3 )* (180/pi), t , y ( : , 3 ) * (180/pi), '--')
xlabel('Time (s)'); ylabel('$\theta$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 4 ) ;
plot ( t , y2 ( : , 4 )* (180/pi),  t , y ( : , 4 ) * (180/pi), '--')
xlabel('Time (s)'); ylabel('$qc/V$ (deg)','Interpreter','latex')

subplot ( 5 , 1 , 5 ) ;
nz2 = (V / g) * ((V/c)*y2(:, 4)' - A(2, :) * y2');
nz1 = (V / g) * ((V/c)*y(:, 4)' - A(2, :) * y');
plot(t, nz2, t, nz1, '--');
xlabel('Time (s)'); ylabel('$n_z$ (-)','Interpreter','latex')

