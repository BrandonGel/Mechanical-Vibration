clc
clear all;
close all;
%%
N = 100;
fs = 100;
dt = 1/fs;


t = (0:N-1).'*dt;
%x = rand([N,1]);
x = 11*sin(2*pi*11*t);
figure(1);
plot(t,x,'linewidth',3);
xlabel('Time [s]');
ylabel('Amplitude [V]');

X = fft(x)*dt;
%%
df = fs/N;
f = (0:N-1).'*df;

figure(2);
subplot(2,1,1);
plot(f,abs(X));
xlabel('Frequency [Hz]');
ylabel('Amplitued [V/Hz]');

subplot(2,1,2);
plot(f,angle(X));
xlabel('Frequency [Hz]');
ylabel('Phase [rad]');
%%
[Gxx,f_half,~,~] = spectral_average(x,fs,1);
figure(3);
plot(f_half,Gxx);