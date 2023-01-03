%Team B2 Hw4
%3/2/21

clc
clear all
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%-Question 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
m = 0.5; % Mass in [kg]
C = 1.2665e-6; %Spring Compliance in [m/N]
R = 10;% Damping Coefficent in [Ns/m]
freqL  = 10;
freqH  = 10000;

%Y(jw) - Mechanical Admittance in Laplace Domain
w = (freqL:freqH);
ss = (1j)*w;
Y = 1./(m.*ss+1./(C.*ss)+R);

%Y(s) plot
mag = 20*log10(abs(Y));
phase = angle(Y)*180/pi;

figure(1)
set(gcf,'position',[0,0,1080,1080])

subplot(2,1,1);
mm = semilogx(w,mag);
set(mm, 'LineWidth', 3)
xlim([freqL,freqH])
ylim([-1e2,0])
grid on
xlabel('Frequency [Hz]', 'FontSize',16,'FontWeight','bold')
ylabel('Magnitude [dB]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('|Y(s)|','FontSize',14);
title('Magnitude Plot of Admittance Y(s)', 'FontSize', 16, 'FontWeight','bold')

subplot(2,1,2);
pp = semilogx(w,phase);
set(pp, 'LineWidth', 3)
xlim([freqL,freqH])
ylim([-1e2,1e2])
grid on
xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Phase [deg]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('\angleY(s)', 'FontSize',14,'Location','best')
title('Phase Plot of  Admittance Y(s)', 'FontSize', 16, 'FontWeight', 'bold')
set(gca,'YTick',(-2:2)*45)
saveas(gcf,'Y.png')
%%
%Define s
s = tf('s');

Y = 1/(m*s+1/(C*s)+R);
wn = sqrt(1/(C*m));
coeff = [25, 100, 20, 20, 20];
wf = [1, 10*wn, 2*wn, wn+1,wn];
xl = [0.05,0.05,0.05,0.05,0.05];
 
for ii = 1:length(wf)
    
    %Solving v(t) & x(t)    
    N = round(coeff(ii)*wn); % Numbers of pts
    dt = 1/(coeff(ii)*wn); %differential time step
    t = 5*(0:N-1)*dt;

    F = 2*sin(wf(ii)*t);

    v = lsim(Y,F,t);
    x = lsim(Y/s,F,t);

    %Ploting v(t)
    figure(ii)
    set(gcf,'position',[0,0,1080,360])
    subplot(1,2,1);
    plot(t,v, 'LineWidth', 1);
    grid on
    xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
    ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    xlim([0,xl(ii)])
    legend('w_f =' + string(wf(ii)),'FontSize',14);
    title('Velocity v(t)', 'FontSize', 16, 'FontWeight','bold')

    %Ploting x(t)
    subplot(1,2,2);
    plot(t,x, 'LineWidth', 1);
    grid on
    xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
    ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    xlim([0,xl(ii)])
    legend('w_f =' + string(wf(ii)),'FontSize',14);
    title('Position x(t)', 'FontSize', 16, 'FontWeight','bold')
    saveas(gcf, string(wf(ii))+'.png');
    
end


%% Changing R & C
%Define s
s = tf('s');

%Parameters
m = 0.5; % Mass in [kg]
C = 1.2665e-6; %Spring Compliance in [m/N]
R = 10;% Damping Coefficent in [Ns/m]

%Y(jw) - Mechanical Admittance in Laplace Domain
w = (freqL:freqH);
ss = (1j)*w;

% set-up
len = 5;
Q = (1:len);
mag=zeros(length(Q),length(w));
phase=zeros(length(Q),length(w));
txt = strings(length(Q),1);
RR = zeros(length(Q),1);
CC = zeros(length(Q),1);
base = 10;

%DSP set-up
wn = sqrt(1/(C*m));
N = round(2.5*wn)*20; % Numbers of pts
dt = 1/(2.5*wn*20); %differential time step
t = (0:N-1)*dt;
wf = 2*wn; %rad/s
F = 2*sin(wf*t);

% x(t) and v(t) arrays
v = zeros(len, N);
x = zeros(len, N);

V=zeros(length(Q),1);
X=zeros(length(Q),1);

%Solving
for ii = 1:length(Q)
    %Y(s) plot
    RR(ii) = R*base^(Q(ii)-3);
    Y = 1./(m.*ss+1./(C.*ss)+RR(ii));
    YY = 1/(m*s+1/(C*s)+RR(ii));
    mag(ii,:) = 20*log10(abs(Y));
    phase(ii,:) = angle(Y)*180/pi;
    txt(ii) = strcat(string(base^(Q(ii)-3)));
    
    %Solving v(t) & x(t)
    v(ii,:) = lsim(YY,F,t);
    x(ii,:) = lsim(YY/s,F,t);
end

%---R Mag & Angle Plot---%
figure(4)
set(gcf,'position',[0,0,1080,1080])

    %----R Mag Plot----%
    subplot(2,1,1);
    mm = semilogx(w,mag);
    set(mm, 'LineWidth', 2)
    xlim([freqL,freqH])
    ylim([-1e2,0])
    grid on
    xlabel('Frequency [Hz]', 'FontSize',16,'FontWeight','bold')
    ylabel('Magnitude [dB]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    leg = legend(txt,'FontSize',14);
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio R_i/R')
    title('Magnitude Plot of Admittance Y(s) with Varying R', 'FontSize', 16, 'FontWeight','bold')

    %----R Angle Plot----%
    subplot(2,1,2);
    pp = semilogx(w,phase);
    set(pp, 'LineWidth', 2)
    xlim([freqL,freqH])
    ylim([-1e2,1e2])
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('Phase [deg]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    leg = legend(txt, 'FontSize',14,'Location','best');
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio R_i/R')
    title('Phase Plot of  Admittance Y(s) with Varying R', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'YTick',(-2:2)*45)
    saveas(gcf,'R.png')
%%
% %----R Velocity Plot----%
figure(5)
set(gcf,'position',[0,0,720,360])
plot(t,v, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,0.01])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio R_i/R')
title('Velocity v(t) with varying R at \omega_f = 1 rad/s', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Rv.png');
%%
%----R Position Plot----%
figure(6)
set(gcf,'position',[0,0,720,360])
plot(t,x, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
xlim([0,0.01])
htitle = get(leg,'Title');
set(htitle, 'String','Ratio R_i/R')
title('Position x(t) with varying R at \omega_f = 1 rad/s', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Rx.png');


%%
%--------------Varying C---------------------------------%

%DSP set-up
wn = sqrt(1/(C*m));
N = round(2.5*wn)*25^2; % Numbers of pts
dt = 1/(2.5*wn*25); %differential time step
t = (0:N-1)*dt;
base = 100;
wf = 2*wn;
F = 2*sin(wf*t);

v = zeros(len, N);
x = zeros(len, N);

%Solving
for ii = 1:length(Q)
    CC(ii) = C*base^(ii-2);
    Y = 1./(m.*ss+1./(CC(ii).*ss)+R);
    YY = 1/(m*s+1/(CC(ii)*s)+R);
    mag(ii,:) = 20*log10(abs(Y));
    phase(ii,:) = angle(Y)*180/pi;
    txt(ii) = strcat(string(base^(Q(ii)-2)));
    
    %Solving v(t) & x(t)
    v(ii,:) = lsim(YY,F,t);
    x(ii,:) = lsim(YY/s,F,t);
    
end

%---C Mag & Angle Plot---%
figure(7)
set(gcf,'position',[0,0,1080,1080])

    %----C Mag Plot----%
    subplot(2,1,1);
    mm = semilogx(w,mag);
    set(mm, 'LineWidth', 2)
    xlim([freqL,freqH])
    ylim([-1e2,0])
    grid on
    xlabel('Frequency [Hz]', 'FontSize',16,'FontWeight','bold')
    ylabel('Magnitude [dB]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    leg = legend(txt,'FontSize',14, 'Location', 'best');
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio C_i/C')
    title('Magnitude Plot of Admittance Y(s) with Varying C', 'FontSize', 16, 'FontWeight','bold')

    %----C Angle Plot----%
    subplot(2,1,2);
    pp = semilogx(w,phase);
    set(pp, 'LineWidth', 2)
    xlim([freqL,freqH])
    ylim([-1e2,1e2])
    grid on
    xlabel('Frequency [Hz]', 'FontSize', 16, 'FontWeight', 'bold')
    ylabel('Phase [deg]', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
    leg = legend(txt, 'FontSize',14,'Location','best');
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio C_i/C')
    title('Phase Plot of  Admittance Y(s) with Varying C', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'YTick',(-2:2)*45)
    saveas(gcf,'C.png')

%----C Velocity Plot----%
figure(8)
plot(t,v, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,0.025])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio C_i/C')
title('Velocity v(t) with varying C', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Cv.png');

%----C Position Plot----%
figure(9)
plot(t,x, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,0.025])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio C_i/C')
title('Position x(t) with varying C', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Cx.png');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%-Question 2%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 1; %kg
g = 9.81; %N/kg
C = 0.1; %m/N
R = 10; %kg/s
theta = [0,pi/6,pi/4,2*pi/3,pi/2];
wn = sqrt(1/m*C);

s = tf('s');
Y = 1/(m*s+1/(C*s)+R);
wf = 10;
dt = 1/(250*wn);
t = (0:dt:10);

figure(1)
for ii = 1:length(theta)
    F = 10*s/(s^2+wf^2) - m*g*sin(theta(ii))/s;
    X = F*Y/s;
    x = impulse(X,t);
    plot(t,x, 'LineWidth', 2)
    hold on
end
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,10])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('0', '\pi/6','\pi/4','\pi/3','\pi/2','FontSize',14);
title('Position vs Time with Varying \theta', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Q2x.png');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%-Question 3%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s = tf('s');

m = 50; %kg
L = 1; %m
d = 0.05; %m
r = 0.2; %m
G = 8.3*10^10; %Pa
zeta = 0.01;
w_f = 215; %rad/s

J = pi/32*d^4;
k = G*J/L;
C = 1/k;

I = 1/2*m*r^2;
wn = sqrt(1/(I*C));
R = 2*I*wn*zeta;

Y = 1/(I*s+R+1/(C*s));
M = w_f^2/(s^2+w_f^2);

Omega = Y*M;
Theta = Omega/s;

t = (0:0.001:2.5);
omega = impulse(Omega,t);
theta = impulse(Theta,t);

figure()
set(gcf,'position',[0,0,960,540])
subplot(2,1,1);
plot(t,omega, 'LineWidth', 2)
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('50kg','FontSize',14);
title('Angular Velocity vs Time', 'FontSize', 16, 'FontWeight','bold')

subplot(2,1,2);
plot(t,theta, 'LineWidth', 2)
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('50kg','FontSize',14);
title('Angular Displacement vs Time', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'Q3.png');





    