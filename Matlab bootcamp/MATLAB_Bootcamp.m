clc;
clear all;
close all;

% Make a comment

%% Vector and Matrix Manipulation
% x = [10+(1j)*10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
% x = x.';

% make a time vector that starts at 0 and has 20 samples in it
N = 100;
dt = 0.01;
time = (0:N-1).'*dt;

freq = 5;
freq2 = 6;
x = sin(2*pi*freq*time);
y = sin(2*pi*freq2*time);

result = x.*y;

% figure(1);
% plot(time,result);

M = [1,2,3;
    4,5,6;
    7,8,9];

%% System of Equations
X = [1, 1, 0;
    1, -1, -1;
    0, 1, -2];

A = [145;0;0];

V = X\A;

%% For Loops
Q = (0:10000000);

tic;
% x = Q.^3+3*Q;
x = zeros(length(Q),1);
for ii = 1:length(Q)
    x(ii) =  Q(ii).^3+3*Q(ii);
end
toc;
%% Plotting
hFig = figure(1);
hPlot = plot(x,Q,'LineWidth',3,'Color',[0 175 0]/255);
xlabel('X Axis Label (Units)','FontSize',16,'FontWeight','bold');
ylabel('Y Axis Label (Units)','FontSize',16,'FontWeight','bold');
title('Title','FontSize',16,'FontWeight','bold');
xlim([0 10*10^20]);
ylim([0 12*10^6]);
grid on;
set(gca, 'FontSize',14,'GridAlpha',0.5);
legend('Legend Text','FontSize',14);
saveas(gcf,'figure_name.tiff');

%% Transfer Functions
s = tf('s');
new_tf = s/(s^2+(5*2*pi)^2);

N = 10000;
dt = 0.001;
time = (0:N-1).'*dt;
x = impulse(new_tf,time);

% figure(1);
% plot(time,x);

R = 1;
C = 0.001;
m = 1;
F = [1;0];
Z = [m*s+R+1/s/C, -1/s/C;
    -1/s/C, m*s+1/s/C];
V = Z\F;
v = impulse(V,time);

figure(2);
plot(time,v);




