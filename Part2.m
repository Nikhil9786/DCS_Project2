close all; clear all; clc;
%% declaring intial values
frequency = 20;
T = 1/frequency;
Tfinal = 1;
t = (0:T:Tfinal)';
N = length(t);
r = 2*pi;
A = 0.2299;
B = 1;
C = 1.4432;
Nvec = [(A-eye(size(A))) B;C 0]\[zeros(size(A));1];
N_x = Nvec(1);
%% reference state from reference signal
reference_state = N_x*r;
%% new augmentation matrices
A_AUG = [A zeros(size(A));-C 1];
B_AUG = [B;0];
%% K matrix
Ts = 0.6;
Mp = 0.2;
zeta = sqrt((log(Mp))^2/(pi^2+(log(Mp))^2));
zeta_omegaN = 4/Ts;
omega_N = zeta_omegaN/zeta;
omega_D = omega_N*sqrt(1-zeta^2);
Roots_cont = [-zeta_omegaN+omega_D*1i -zeta_omegaN-omega_D*1i];
Roots_disc = exp(Roots_cont*T);
z_star = [Roots_disc(1) Roots_disc(2)];
K_matrix = acker(A_AUG,B_AUG,[z_star]);
%% L matrix
Ts_estimate = 0.2;
zeta_omegaN_estimate = 4/Ts_estimate;
Roots_cont_estimate = [-zeta_omegaN_estimate];
Roots_disc_estimate = exp(Roots_cont_estimate*T);
beta_star = [Roots_disc_estimate];
L_matrix = acker(A', C',beta_star)';

%% Simulation
x = zeros(1,N);
x_hat = zeros(1,N);
y = zeros(1,N);
u = zeros(1,N);
e = zeros(1,N);
xI = zeros(1,N);
error = zeros(1,N);
x_hat(:,1) = zeros(1,1);
x(:,1) = zeros(1,1);
y(1) = C*x(:,1);
e(1) = r-y(1);
xI(1) = 0;
u(1) = 0;
for k=2:N
    x(k) = A*x(k-1)+B*u(k-1);
    y(k) = C*x(k);
    x_hat(k) = (A-L_matrix*C)*x_hat(k-1)+B*u(k-1)+L_matrix*y(k-1);
    e(k) = r-y(k);
    xI(k) = xI(k-1)+e(k-1); 
    u(k) = K_matrix*([reference_state; 0]-[x_hat(k); xI(k)]);
    error(k) = x_hat(k) - x(k);
end

%% Plots
figure(1)
plot(t,y(1,:),'k-o', t,r,'r-o', 'linewi',2, 'MArkersize', 6)
grid on
hold on
title('System Output, y')
xlabel('t')
ylabel('Angular Velocity')
legend('Motor Angular Velocity', 'reference line')

figure(2)
stairs(t,u, 'k-', 'linewi',2)
grid on 
hold on
title('Control Input')
xlabel('t')
ylabel('Voltage')
legend('Input Voltage,[V]')

figure(3)
plot(t,error,'k-','linewi',2)
ylim([-1 1])
grid on
title('State estimation error')
xlabel('t')
ylabel('Amplitude')
legend('Error')
grid on 