clc;
close all;
clear all;

%% Vibration Matlab Topic 3 I

s = tf('s');
t = (0:0.001:10);

m = 0.1;
c = 1/12; 
wn = sqrt(1/(m*c));
wf = [1, 10*wn, 2*wn, wn+1, wn];
x0 = 0;
v0 = 0;

for i = 1:length(wf)
    Z = (m*s+1/(s*c));
    Y = 1/Z;
    Fv0 = m*v0;
    Fx0 = -x0/c/s;
    Fi = s/(s^2+wf(i)^2);
    F = Fi+Fx0 + Fv0;

    V = Y*F;
    X = V/s;
    v = impulse(V,t);
    x = impulse(X,t);
    f = impulse(F,t);
    figure(1)
    plot(t,v, 'LineWidth', 3);
    hold on
    figure(2)
    plot(t,x, 'LineWidth', 3);
    hold on
    txt = int2str(wf(i));
    %legend(txt, 'f')
    
end
%% 3.1.2 Security Camera Mount
s = tf('s');


m = 3; %kg
L = 0.2; %m
E = 200*10^9; %Pa
a = 0.01; %m

I = a^4/12;
k = 3*E*I/L^3;
c = 1/k;
wn = sqrt(1/(m*c));

Z = (m*s+1/(s*c));
Y = 1/Z;

freq = 10;
wf = 2*pi*freq;
f0 = 15;
F = f0*s/(s^2+wf^2);

V = Y*F;
X = V/s;

t = (0:1/(2.5*wn):1);
f = impulse(F,t);
v = impulse(V,t);
x = impulse(X,t);

figure(1)
plot(t,x, 'LineWidth', 3)
legend('X')

fprintf('Natural Frequency: %f Hz\n', wn)
fprintf('Forced Frequency: %f Hz\n', freq)

%% Vibration Matlab Topic 3 II

om = (0.1:0.01:10);
ss = om*1j;

R = 100;
m = 10;
k = 70;
Y = (R+k./ss)./(R+m*ss+k./ss);

figure()
plot(om,abs(Y));


%%
s = tf('s');
vx = [1,2,3];
lamb = 6;
wf = 2*pi*vx/lamb;
YY = (R+k/s)/(R+m*s+k/s);
t = (0:0.1:10);

% for ii = 1:length(wf)
%     X_b = 0.02*wf(ii)/(s^2+wf(ii)^2);
%     X = YY*X_b;
%     x(ii,:) = impulse(X,t);
% end

X_b = 0.02.*wf./(s^2+wf.^2);
X = YY*X_b;
x(ii,:) = impulse(X,t);

figure(2)
plot(t,x, 'LineWidth', 3);
grid on;
xlabel('Time [sec]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Displacement [m]', 'Fontsize', 12, 'FontWeight', 'bold');
legend(string(vx)+'=[m/s]', 'FontSize', 8);
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3, 'Ytick', (-0.02:0.01:0.02));
title('Topic 3 Displacement', 'FontSize', 14', 'FontWeight', 'bold');








