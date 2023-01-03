clc;
clear all;
close all;

%% Question 3: Analysis of Mechanical System

%Parameter
R = 4.35; % in [kg/s]
C1 = 1/1.2; % in [m/N]
C2 = 1/4; % in [m/N]
C3 = 1/0.33; % in [m/N]
m = 0.65; % in [kg]

%Frequency Data Points
omega = 0:0.01:1000; % Hz
ss = 1j*omega;

%Solving for Admittance T.F Y(s)
Ceq = C2+C3;
Z = 1./(Ceq*ss) + 1./(C1*ss) + R + m*ss;
Y = 1./Z;

% Y(s) Plot
figure(1)
semilogx(omega,abs(Y));
xlabel('Angular Frequency (rad/s)');
ylabel('Admittance Magnitude(m/(s*N))');

% Finding the natural frequency and damping ratio
[~,wn] = max(real(Y));
[~,w1] = max(imag(Y));
[~,w2] = min(imag(Y));

omega_n = omega(wn);
omega_1 = omega(w1);
omega_2 = omega(w2);

Q = omega_n/(omega_2-omega_1);
zeta = 1/2/Q;

fprintf('wn = %f rad/s\n',omega_n);
fprintf('zeta = %f rad/s\n',zeta);

%% Question 4: Skyscraper Damping

%Parameter
L = 12.6; % Length of the pendulum in [m]
m = 6.6e5; % Mass in [kg]
C = 8.0e5; % Damping coefficient in [kg/s]
g = 9.81; % gravity in [N/kg]
w_des = 0.620;

%Frequency Data Points
omega = (0:0.001:1000);
ss = 1j*omega;

%Solving for Admittance T.F Y(s)
Z = m*L^2*ss+m*g*L./ss+C*L^2;
Y = 1./Z;

mag = 20*log10(abs(Y));
phase = angle(Y)*180/pi;

figure();
ax1 = subplot(2,1,1);
semilogx(omega,mag);
grid on;
xlabel('Frequency [Hz]', 'Fontsize',16, 'Fontweight', 'bold');
ylabel('Magnitude [m/s/N]', 'Fontsize', 16, 'Fontweight', 'bold');
title('Admittance Y(s) Bode Plot', 'Fontsize', 16, 'Fontweight', 'bold');
set(gca, 'Fontsize', 14, 'GridAlpha', 0.5, 'MinorGridAlpha', 0.1, 'MinorGridLineStyle', '-');
legend('|Y(s)|', 'Fontsize', 14)

ax2 = subplot(2,1,2);
semilogx(omega,phase);
grid on;
xlabel('Frequency [Hz]', 'Fontsize', 16, 'Fontweight', 'bold');
ylabel('Phase [deg]', 'Fontsize', 16, 'Fontweight', 'bold');
set(gca, 'Fontsize', 14, 'GridAlpha', 0.5, 'MinorGridAlpha', 0.1, 'MinorGridLineStyle', '-', 'YTick', (-2:2)*45);
legend('\angle Y(s)', 'Fontsize', 14)
saveas(gcf,'Q4Plot.png')

[~,omega_n] = max(mag);
wn = omega(omega_n); %Natural Frequency in [rad/s]

% Finding the natural frequency and damping ratio
[~,wn] = max(real(Y));
[~,w1] = max(imag(Y));
[~,w2] = min(imag(Y));

omega_n = omega(wn);
omega_1 = omega(w1);
omega_2 = omega(w2);

Q = omega_n/(omega_2-omega_1);
zeta = 1/2/Q;

fprintf('wd = %f rad/s\n',omega_n*sqrt(1-zeta^2));
