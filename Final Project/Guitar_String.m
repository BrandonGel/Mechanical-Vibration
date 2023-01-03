
clc;
clear all;
close all;
s = tf('s');

%Constants
%String Guage[m]:
r = [0.2540, 0.3302, 0.4310, 0.6604, 0.9144, 1.1684]*10^-3; %E(plain), B(plain), G(plain), D(Wound), A(wound), E(wound)
E = 200E9;% Steel Young Modulus
rr=r(1);
rho = 8000; %Steel Density
R = 0;

%Number of nodes:
%n = input("Type in an odd numbers of interior nodes: ");
n = 15;

%Type in the length of string:
%l = input("\nType in the length of the string: ");
l = 10;

%Derived Parameter
dl = l/n;
A = pi*rr^2;
k = E*A/dl;
R = R/n;
m = rho*A*dl;

%Location of the Interior and Boundary Nodes
loc = zeros(n+2,1);
loc(1) = 0;
loc(n+2) = l;
for ii = 2:n+1
    loc(ii) = (ii-1.5)*dl;
end

%Z,Y,H T.F
%Interiors Nodes are evenly spaced
%Distance b/w boundary and interior nodes are 1/2 of the interior nodes'

if n == 1
    Z = m*s+k/s+R;
else
    %For 1st boundary node
    Z(1,1) = m*s + (1.5*k)/s + 1.5*R;
    Z(1,2) = -(k/s + R);
    for ii = 2:n-1
        Z(ii,ii-1) = Z(ii-1,ii);
        Z(ii,ii) = m*s -2*Z(ii,ii-1);
        Z(ii,ii+1) = Z(ii,ii-1);
    end   
    %For last boundary node
    Z(n,n-1) = Z(1,2);
    Z(n,n) = Z(1,1); 
end
Z = minreal(Z);
Z = zpk(Z);
Y = inv(Z);
Y = minreal(Y);
H = Y/s;

%Pluck length at the center
xpluck = .05;
x0 = zeros(n,1);
for ii = 1:(n+1)/2
   x0(ii) = 2*xpluck/l*loc(ii+1);
   x0(n+1-ii) = x0(ii);
end

%Spring Force Calculation
f = zeros(n+1,1);
f(1) = k/2*x0(1);
f(n+1) = -k/2*x0(n);
for ii = 2:n
    f(ii) = k*(x0(ii)-x0(ii-1));
end

%Initial Force Voltage
for ii = 1:n
    F(ii,1) = (f(ii+1)-f(ii))/s;
end

%%
%Postprocessing
leg = strings(n);
for ii = 1:n
   leg(ii) = "x" + int2str(ii); 
end

%Calculating the transient response
wn = abs(max(pole(Y)));
X = H*F+x0/s;
X =minreal(X);
t = (0:1/(2.5*wn):0.1);
x = impulse(t,X);

%%
%Plotting the transient response
figure(1);
plot(t,x,'linewidth', 3)
ylabel('displacment [m]', 'fontsize', 12, 'fontweight', 'bold');
xlabel('Time [s]', 'fontweight', 'bold')
legend(leg);

len = int16(length(t)/2);
pos = zeros(len,n+2);
for ii = 1:len
    figure(2)
    pos(ii,2:n+1) = x(ii,:);
    plot(loc,pos(ii,:),'linewidth',3);
    ylim([-xpluck xpluck])
    %pause(0.05);
    %hold on;
    %val = input('Press Enter to continue to next Iteration');
end


%%
om = (100:1:10000);
[mag,phase,om]=bode(H,om);
figure(2)
for ii = 1:n
    for jj = 1:n
        subplot(n,n,n*(ii-1)+jj)
        semilogx(om,mag(ii,jj),'linewidth', 3')
        %str = 'x' + strings(ii) + strings(jj);
        ylabel('disp', 'fontsize', 12, 'fontweight', 'bold');

    end
end
