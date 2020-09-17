clear all; close all; clc;
%% Declaring intial values
Tfinal = 1.5;
frequency = 100;
T = 1/frequency;
t = (0:T:Tfinal)';
N = length(t);
n = 5768280;
alpha = 115856;
gamma = 19230;
Ts = 1.5;
Mp = 0.3;
zeta = sqrt((log(Mp))^2/(pi^2+(log(Mp))^2));
zeta_omegaN = 4/Ts;
omega_N = zeta_omegaN/zeta;
omega_D = omega_N*sqrt(1-zeta^2);
Roots_cont = [-zeta_omegaN+omega_D*1i -zeta_omegaN-omega_D*1i];
Roots_disc = exp(Roots_cont*T);
z_star = [Roots_disc(1) Roots_disc(2)];
Ac = [0 1;n/alpha 0];   %continous state space model
Bc = [0 -gamma/alpha]';
A = expm(Ac*T);         %continous to discrete state space
B = (A - eye(length(A)))*Ac^-1*Bc;
A(1:2,3) = 0;
A(3,1:2) = 0;
A(3,3) = 1;
B(3,1) = T;
K_matrix = acker(A,B,[Roots_disc(1), Roots_disc(2), abs(Roots_disc(1))]); %calculating k matrix
C = [0 1 0;0 0 1];
Ts_estimate = 0.5;
Mp_estimate = 0.2;
zeta_estimate = sqrt((log(Mp_estimate))^2/(pi^2+(log(Mp_estimate))^2));
zeta_omegaN_estimate = 4/Ts_estimate;
omega_N_estimate = zeta_omegaN_estimate/zeta_estimate;
omega_D_estimate = omega_N_estimate*sqrt(1-zeta_estimate^2);
Roots_cont_estimate = [-zeta_omegaN_estimate+omega_D_estimate*1i ...
    -zeta_omegaN_estimate-omega_D_estimate*1i];
Roots_disc_estimated = exp(Roots_cont_estimate*T);
beta_star = [Roots_disc_estimated(1) Roots_disc_estimated(2)];
char_Equation = poly([beta_star]);
L_matrix = place(A', C',[Roots_disc_estimated(1) Roots_disc_estimated(2) ...
    abs(Roots_disc_estimated(1))])';
%% Simulation
x = zeros(3,N);
x_hat = zeros(3,N);
y = zeros(2,N);
y_hat = zeros(2,N);
u = zeros(1,N);
x(:,1) = [0;0.2;0];
x_hat(:,1) = [0;0;0];
y_hat(1) = C(1,:)*x_hat(:,1);
u(1) = -K_matrix*x(:,1);
u(1) = 0;   
error(:,1) = x_hat(:,1)-x(:,1);
for k = 2:N
    x(:,k) = A*x(:,k-1)+B*u(k-1);
    x_hat(:,k) = A*x_hat(:,k-1) + B*u(k-1) - L_matrix * (y_hat(:,k-1)-y(:,k-1));
    y_hat(:,k) = C(:,:)*x_hat(:,k);
    y(:,k) = C(:,:)*x(:,k);
    u(k) = -K_matrix*x_hat(:,k);
    error(:,k) = x_hat(:,k) - x(:,k);
end

%% Plots
figure(1)
plot(t,y(1,:)*180/pi, 'k-o', t,y(2,:)*180/pi,'r-o', 'linewi',1,'MarkerSize',4)
grid on
title('Measured Angle')
xlabel('Time')
ylabel('\theta/s')
legend('pitch rate','wheel angular veliocity')

figure(2)
stairs(t,u*180/pi,'k-','linewi',1, 'MArkersize',4)
grid on
title('Angular Acceleration')
xlabel('Time')
ylabel('\theta/s^2')
legend('wheel angular acceleration')

figure(3)
plot(t,error(1,:),'k-*', t,error(2,:),'r-*', t,error(3,:),'g-*','linewi',1,'MarkerSize',4)
grid on
title('State Estimation Error')
xlabel('Time')
ylabel('Error')
legend('Error 1','Error 2','Error 3')