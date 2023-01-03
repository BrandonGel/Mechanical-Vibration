%Team B2 Hw3
%2/13/21

clc
clear all
close all

%%
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

figure()
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

% IC's
x0 = 0.3 ; %Inital Displacement of the Mass [m]
v0 = 1; %Inital Velocity of the Mass [m/s]

%V(s) IC's in Laplace Domain
FX0 = -1/C*x0/s;
FV0 = m*v0;

%X(s) IC's in Laplace Domain
X0 = x0/s;

%TF: V(s) = F(s)*Y(s)
F = FX0+FV0;
Y = 1/(m*s+1/(C*s)+R);
V = F*Y;

%TF: X(s)-X(0) = integrate(V(s)) = V(s)/s
X = V/s+X0;

%Solving v(t) & x(t)
wn = sqrt(1/(C*m));
N = round(2.5*wn)/2; % Numbers of pts
dt = 1/(2.5*wn); %differential time step
t = (0:N-1)*dt;
v = impulse(V,t);
x = impulse(X,t);

%Ploting v(t)
figure()
plot(t,v, 'LineWidth', 1, 'color', 'r');
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('v(t)','FontSize',14);
title('Velocity v(t)', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'v.png');

%Ploting x(t)
figure()
plot(t,x, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('x(t)','FontSize',14);
title('Position x(t)', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'x.png');

%% Changing R & C
%Define s
s = tf('s');

%Parameters
m = 0.5; % Mass in [kg]
C = 1.2665e-6; %Spring Compliance in [m/N]
R = 10;% Damping Coefficent in [Ns/m]

% IC's
x0 = 0.3 ; %Inital Displacement of the Mass [m]
v0 = 1; %Inital Velocity of the Mass [m/s]

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
base = 5;

%V(s) IC's in Laplace Domain
FX0 = -1/C*x0/s;
FV0 = m*v0;
F = FX0+FV0;

%X(s) IC's in Laplace Domain
X0 = x0/s;

% x(t) and v(t) arrays
v = zeros(len, N);
x = zeros(len, N);
%---------------------------------------------------------%

%DSP set-up
wn = sqrt(1/(C*m));
N = round(2.5*wn)/2; % Numbers of pts
dt = 1/(2.5*wn); %differential time step
t = (0:N-1)*dt;

V=zeros(length(Q),1);
X=zeros(length(Q),1);

%Solving
for ii = 1:length(Q)
    %Y(s) plot
    RR(ii) = R*base^(Q(ii)-3);
    Y = 1./(m.*ss+1./(C.*ss)+RR(ii));
    mag(ii,:) = 20*log10(abs(Y));
    phase(ii,:) = angle(Y)*180/pi;
    txt(ii) = strcat(string(base^(Q(ii)-3)));
    
    %TF: V(s) = F(s)*Y(s)
    %TF: X(s)-X(0) = integrate(V(s)) = V(s)/s
    V = F/(m*s+1/(C*s)+RR(ii));
    X = V/s+X0;
    
    %Solving v(t) & x(t)
    v(ii,:) = impulse(V,t);
    x(ii,:) = impulse(X,t);
end

%---R Mag & Angle Plot---%
figure()
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
    leg = legend(txt, 'FontSize',14,'Location','best')
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio R_i/R')
    title('Phase Plot of  Admittance Y(s) with Varying R', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'YTick',(-2:2)*45)
    saveas(gcf,'R.png')

%----R Velocity Plot----%
figure()
plot(t,v, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio R_i/R')
title('Velocity v(t) with varying R', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'v2.png');

%----R Position Plot----%
figure()
plot(t,x, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio R_i/R')
title('Position x(t) with varying R', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'x2.png');





%--------------Varying C---------------------------------%

%DSP set-up
wn = sqrt(1/(C*m));
N = round(2.5*wn)/2*25; % Numbers of pts
dt = 1/(2.5*wn*25); %differential time step
t = (0:N-1)*dt;
base = 100;

v = zeros(len, N);
x = zeros(len, N);

%Solving
for ii = 1:length(Q)
    CC(ii) = C*base^(ii-2);
    Y = 1./(m.*ss+1./(CC(ii).*ss)+R);
    mag(ii,:) = 20*log10(abs(Y));
    phase(ii,:) = angle(Y)*180/pi;
    txt(ii) = strcat(string(base^(Q(ii)-2)));
    
    %V(s) IC's in Laplace Domain
    FX0 = -1/CC(ii)*x0/s;
    F = FX0+FV0;
    
    %TF: V(s) = F(s)*Y(s)
    V = F/(m*s+1/(CC(ii)*s)+R);
    X = V/s+X0;
    
    %Solving v(t) & x(t)
    v(ii,:) = impulse(V,t);
    x(ii,:) = impulse(X,t);
    
end

%---C Mag & Angle Plot---%
figure()
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
    leg = legend(txt, 'FontSize',14,'Location','best')
    htitle = get(leg,'Title');
    set(htitle, 'String','Ratio C_i/C')
    title('Phase Plot of  Admittance Y(s) with Varying C', 'FontSize', 16, 'FontWeight', 'bold')
    set(gca,'YTick',(-2:2)*45)
    saveas(gcf,'C.png')

%----C Velocity Plot----%
figure()
plot(t,v, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Velocity [m/s]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,0.5])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio C_i/C')
title('Velocity v(t) with varying C', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'v3.png');

%----C Position Plot----%
figure()
plot(t,x, 'LineWidth', 1);
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Position [m]', 'FontSize', 16, 'FontWeight', 'bold')
xlim([0,0.5])
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
leg = legend(txt,'FontSize',14);
htitle = get(leg,'Title');
set(htitle, 'String','Ratio C_i/C')
title('Position x(t) with varying C', 'FontSize', 16, 'FontWeight','bold')
saveas(gcf,'x3.png');


    