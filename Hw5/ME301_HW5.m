%Team B2 Hw5
%3/6/21

clc
clear all
close all

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%-Question 1%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
M = 81; % Washing Machine Mass in [kg]
m = 16; % Cloth Mass in [kg]
k = 20e4; % Spring Constant in [N/m]
R = 10;% Damping Coefficent in [Ns/m]
r = 0.3;% Radius of the Washing Machine in [m]
wf = 900 *2*pi/60; % Angular velocity produced by the Washing Machine in [rad/s]

C = 1/k; % Spring Compliance in [N/m]
wn = sqrt(1/((m+M)*C)); % Natural Frequency of the System in [rad/s]

%Define s
s = tf('s');


%TF of the Washing Machine
Z = [1, 0; -m*s, (M+m)*s+R+1/(C*s)];
%Velocity Input
V = [r*wf*s/(s^2+wf^2); 0];

V_1 = Z\V;
X_1 = V_1/s;
t = (0:1/(wn*10):60).';
x_1 = 100*impulse(t,X_1); %Washing Machine Displacment in [cm]


%Y(s) plot
figure(1)
set(gcf,'position',[0,0,1080,720])
plt = plot(t,x_1(:,2),'LineWidth', 3);
set(plt, 'LineWidth', 3)
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Displacment [cm]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend('y(t)','FontSize',14);
title('Washing Machine Displacement Plot', 'FontSize', 16, 'FontWeight','bold')
print(gcf,'Disp1.png', '-dpng', '-r300');

xx = x_1(round(50*(wn*10)):end,2);
fprintf('Max steady-state displacement: %fcm', max(xx));
fprintf(newline);
%%
Rg = 10; % Ground Damper in [kg/s]
kg = 100; % Ground Spring Constant in [N/m]
Cg = 1/kg; % Ground Spring Compliance in [m/N]
mg = 19.5;
txt = "R_g = " + string(Rg) +'kg/s'+ newline  +"k_g = " + string(kg) +'N/m' + newline  +"m_g = " + string(mg) + 'kg';
%C2 = 1/(wf^2*(m+M))-C;
%TF of the Washing Machine
Z = [1, 0, 0;
    -m*s, (M+m)*s+R+1/(C*s), -(R+1/(C*s));
    0, -(R+1/(C*s)), (R+1/(C*s)+Rg+1/(Cg*s)+mg*s)];

V = [r*wf*s/(s^2+wf^2); 0; 0];

V_1 = Z\V;
X_1 = V_1/s;
t = (0:0.001:60).';
x_2 = 100*impulse(t,X_1); %Washing Machine Displacment in [cm]


%Y(s) plot
figure(2)
set(gcf,'position',[0,0,1080,720])
plt = plot(t,x_2(:,2),'LineWidth', 3);
set(plt, 'LineWidth', 3)
grid on
xlabel('Time [s]', 'FontSize',16,'FontWeight','bold')
ylabel('Displacment [cm]', 'FontSize', 16, 'FontWeight', 'bold')
set(gca, 'FontSize',14,'GridAlpha',0.5,'MinorGridAlpha', 0.5);
legend(txt,'FontSize',14);
title('Washing Machine Displacement Plot', 'FontSize', 16, 'FontWeight','bold')
print(gcf,'Disp2.png', '-dpng', '-r300');

xx = x_2((50/0.001:end),2);
fprintf('Max steady-state displacement: %fcm', max(xx));
fprintf(newline);
