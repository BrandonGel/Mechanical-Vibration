clc;
clear all
close all
%%
k = 2e3;
m = 10;
R = 20;

s = tf('s');
N = 10000;
dt = 0.001;
t = (0:N-1).'*dt;
wn = sqrt(k/m);
wf = 10;
t0 = 1;

Y = 1/(R+m*s+k/s);
H = Y/s;

h = impulse(H,t);  %Impulse Response Fucntion
%% F1 & f1

%Laplace Transform
F1 = exp(-s*t0);
V1 = Y*F1;
X1 = V1/s;
f1 = impulse(F1,t);
x1 = impulse(X1,t);

%lsim
% f1 = zeros(N,1);
% f1(t==t0) = 1;
% x1 = lsim(H,f1,t);

figure(1)
subplot(2,1,1);
plot(t,f1);

subplot(2,1,2);
plot(t,x1);

%% F2 & f2

%LaPlace Transfrorm
F2 = wf/(s^2+wf^2);
V2 = Y*F2;
X2 = V2/s;
%f2 = impulse(F2,t);
%x2 = impulse(X2,t);

% %lsim
f2 = zeros(N,1);
f2(t == t0) = 1;
f2 = sin(wf*t) + f2;
x2 = lsim(H,f2,t);

X = H*F2;
figure(1)
subplot(2,1,1);
plot(t,f2);

subplot(2,1,2);
plot(t,x2);
%% F3 & f3

%LaPlace Transfrorm
F3 = 1/s*exp(-t0*s);
% V3 = Y*F3;
% X3 = V3/s;
% f3 = impulse(F3,t);
% x3 = impulse(X3,t);

%lsim
%f3 = u(t-t0) + sin(wf*t);
f3 = zeros(N,1);
f3(t >= 1) = 1;
x3 = lsim(H,f3,t);

figure(1)
subplot(2,1,1);
plot(t,f3);

subplot(2,1,2);
plot(t,x3);
%% F4 & f4
f4 = t;
f4(t >= t0) = sin(wf*(t(t>=t0)));
x4 = lsim(H,f4,t);

figure(1)
subplot(2,1,1);
plot(t,f4);
ylim([-1.5,1.5]);

subplot(2,1,2);
plot(t,x4);
%% Examine X,F, and H in the LaPlace (Frequency Domain)
w = [0.1:0.001:1000];
[Hmag, Hphase] = bode(H,w);
[Fmag, Fphase] = bode(F2,w);
[Xmag, Xphase] = bode(X,w);

Hmag = squeeze(Hmag);
Fmag = squeeze(Fmag);
Xmag = squeeze(Xmag);

Hphase = squeeze(Hphase);
Fphase = squeeze(Fphase);
Xphase = squeeze(Xphase);

figure(101)
subplot(2,1,1);
semilogx(w,20*log10(Hmag), w,20*log10(Fmag),w,20*log10(Xmag))
legend('H','F','X')

subplot(2,1,2);
semilogx(w,Hphase, w,Fphase,w,Xphase)
legend('H','F','X')

%% Examine H, Y, and I in the LaPlace (Frequency Domain)
w = [0.1:0.001:1000];
H = Y/s;
I = s*Y;

[Hmag, Hphase] = bode(H,w);
[Ymag, Yphase] = bode(Y,w);
[Imag, Iphase] = bode(I,w);

Hmag = squeeze(Hmag);
Ymag = squeeze(Ymag);
Imag = squeeze(Imag);

Hphase = squeeze(Hphase);
Yphase = squeeze(Yphase);
Iphase = squeeze(Iphase);

figure(101)
semilogx(w,20*log10(Hmag), w,20*log10(Ymag),w,20*log10(Imag))
legend('H','F','X')


h = impulse(H,t);
y = impulse(Y,t);
i = impulse(I,t);

figure(3)
plot(t,y*10,t,h*100,t,i,'linewidth', 2);
legend('y','h','i')
[~,Hloc] = max(Hmag);
[~,Yloc] = max(Ymag);
[~,Iloc] = max(Imag);

w(Hloc)
w(Yloc)
w(Iloc)



