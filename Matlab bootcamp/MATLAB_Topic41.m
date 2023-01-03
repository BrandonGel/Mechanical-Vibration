
clc
clear all
close all
%%

%Constants & Parameters
m1 = 10;
m2 = 15;
k1 = 20;
k2 = 24;
k3 = 22;
R1 = 0;
R2 = 0;
R3 = 0;

s = tf('s');
Z = [ m1*s+R1+R2+k1/s+k2/s, -(R2+k2/s);
    -(R2+k2/s), m2*s+R2+R3+k2/s+k3/s;];
Y = inv(Z);
H = Y/s;


om = [0.1: 0.001: 100];
figure(1)
bode(H,om)
[mag, phase, w] = bode(H,om);
%%
x01 = 0.01;
x02 = 0.015;
x0 = [x01; x02];

F1 = k1*x01/s;
F2 = k2*(x02-x01)/s;
F3 = k3*x02/s;
F = [-F1+F2; -F2-F3];
X = H*F;
X = X+x0/s;
t = (0:0.001:100);
x = impulse(X,t);

figure(2)
subplot(2,1,1);
plot(t,x(:,1),'linewidth', 3')
ylabel('X_{1}(t) [m]', 'fontsize', 12, 'fontweight', 'bold');
subplot(2,1,2);
plot(t,x(:,2),'linewidth', 3')
ylabel('X_{2}(t) [m]', 'fontsize', 12, 'fontweight', 'bold');
xlabel('Time [s]', 'fontsize', 12, 'fontweight', 'bold')

figure(3)
plot(t,x,'linewidth', 3')
ylabel('X_{2}(t) [m]', 'fontsize', 12, 'fontweight', 'bold');
xlabel('Time [s]', 'fontweight', 'bold')

%%
[reals, imags] = nyquist(H,w);
H_w = reals + 1j*imags;


fig = figure(1);
subplot(2,2,1)
semilogx(w,20*log10(abs(squeeze(H_w(1,1,:)))));

subplots(2,2,2)
semilogx(w,20*log10(abs(squeeze(H_w(1,2,:)))));

subplots(2,2,3)
semilogx(w,20*log10(abs(squeeze(H_w(2,1,:)))));

subplots(2,2,4)
semilogx(w,20*log10(abs(squeeze(H_w(2,2,:)))));

%%
[~,mode_inds] = findpeaks(squeeze(abs(H_w(1,1,:))));
modes = w(mode_inds);

%%
v = imag(squeeze(H_w(1,1:2,:)).');
v = v(mode_inds,:).';

figure(7)
plot(v+[-5,-5;5,5],[1,2;1,2],'LineStyle', 'none', 'Marker', 'o');
xlim([-20,20]);
ylim([0,3]);

%%
H_w11 = squeeze(abs(H_w(1,1,:)));
H_w12 = squeeze(abs(H_w(1,2,:)));
H_w21 = squeeze(abs(H_w(2,1,:)));
H_w22 = squeeze(abs(H_w(2,2,:)));
[pks11,mod_inds11] = findpeaks(H_w11);
[pks12,mod_inds12] = findpeaks(H_w12);
[pks21,mod_inds21] = findpeaks(H_w21);
[pks22,mod_inds22] = findpeaks(H_w22);
modes11 = w(mod_inds11);
modes12 = w(mod_inds12);
modes21 = w(mod_inds21);
modes22 = w(mod_inds22);

figure(4)
plot(w,imags(2,2));





