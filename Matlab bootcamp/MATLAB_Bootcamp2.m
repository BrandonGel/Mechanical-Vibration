clc;
clear all;
close all;

R = 1 % kg/s
C = 1; % m/N
m = 1; %kg

x0 = 0.01; % m
v0 = 0.01; % m/s

%% Time Domain
omega_n = sqrt(1/(m*C));
zeta = R/(2*m*omega_n);
omega_d = omega_n*sqrt(1-zeta^2);

N = 100000;
dt = 0.001;
t = (0:N-1).'*dt;

x = (x0*cos(omega_d*t)+(v0+zeta*omega_n*x0)/omega_d*sin(omega_d*t)).*exp(-zeta*omega_n*t);

figure(1);
plot(t,x);
xlabel('Time [s]');
ylabel('Displacement [m]');

%% LaPlace Domain
s = tf('s');

X = x0*(s+zeta*omega_n)/((s+zeta*omega_n)^2 + omega_d^2) + (v0+zeta*omega_n*x0)/omega_d*omega_d^2/((s+zeta*omega_n)^2+omega_d^2);
x_s = impulse(X,t);

V = s*X-x0;
v_s = impulse(V,t)

figure(2)
plot(t,x,t,x_s);
xlabel('Time [s]');
ylabel('Displacement [m]');

figure(3)
plot(t,v_s);
xlabel('Time [s]');
ylabel('Displacement [m]');


%% LaPlace Domain with Impedance Matrices (Admittance)

Z = R + 1/s/C + s*m;
Y = 1/Z  % 1/(m*s_1/s/C+R)

F = -1/C/s*x0 + m*v0;

V = Y*F;
v_s2 = impulse(V,t);

X = V/s + x0/s;
x_s2 = impulse(X,t);

figure (4)
plot(t, v_s2, t, x_s2);

%% LaPlace Domain -> Frequency Domain
% s = j*omega
% freq = 0:0.01:100; % Hz
% omega = 2*pi*freq; % rad/s
omega = 0:0.0001:100; % Hz

Y_omega = 1./(m*(1j*omega) + 1./(1j*omega)/C + R);

figure(6)
plot(omega,abs(Y_omega));
xlabel('Angular Frequency (rad/s)');
ylabel('Admittance Magnitude(m/(s*N))');

figure(7)
subplot(2,1,1);
semilogx(omega,real(Y_omega));
xlabel('Angular Frequency (rad/s)');
ylabel('Admittance Magnitude(m/(s*N))');

subplot(2,1,2);
semilogx(omega,imag(Y_omega));
xlabel('Angular Frequency (rad/s)');
ylabel('Admittance Imag Magnitude(m/(s*N))');

[~,ind_wn] = max(real(Y_omega));
[~,ind_w1] = max(imag(Y_omega));
[~,ind_w2] = min(imag(Y_omega));

omega_n2 = omega(ind_wn);
omega_1 = omega(ind_w1);
omega_2 = omega(ind_w2);

Q = omega_n/(omega_2-omega_1);
zeta = 1/2/Q;



