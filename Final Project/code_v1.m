clc;
clear all;
close all;

%% ============================= Parameters ===============================

% ------------ Guitar String ---------
% Length
Length = [0.6477, 10];
ll = 2;
strLength = Length(ll); % meters

% Diameter
strDiameter = 0.001; %meters

% Density
strRho = 8.90*1000; %kg/m^3

% Mass of String
strMass = pi*strDiameter^2/4*strLength*strRho; %kg
        % Area * Length * Density


%------------Damping and Compliance---------
% Damping
R = 0; %i forgot the units

% Time matrix stuff
s = tf('s');
N = 100000;
dt = 0.001;
t = (0:N-1).'*dt;
% Omegas from 0.1 to 10 in 0.001 increments
omega = (10:0.1:1000);

% ==================== User Defined/Calculated Parameters ================
% Number of Nodes
N = 5;

% Distance Between Nodes
nodeDist = strLength/(N+1); %m

% Compliance Between Node
Tension = 80; 
k = pi^2* Tension/nodeDist;
C = 1/k; %idk lmao, based on string tension, use SQRT(Tension/LinearDensity) to find SoS


%Tension to Compliance
%More Tension, higher frequency, more stiff, less compliant
%Therefore, c = frequency of sound * wavelength of wave
%              fn = 1/2pi * SQRT(k/M)  wavelength 2L/n, n= the number of the harmonic.
%                                      Lambda= 2L
%              



% Mass of each Node
nodeMass = strMass/N; %kg

% ---------- Impedance Calculations ---------
% Mass Impedance
Zm = nodeMass*s;
 
% Spring Impedance
Zk = 1/(C*s);
 
% Damping Impedance
Zr = R/(N+1);

%----------- DEBUG ONLY --------------
% Zm = 1;
% Zk = 2;
% Zr = 5;

% === oh god impedance save me but lucky for me hannah had a solution ====
% The first and last entries are speshul because they aren't bracketed by
% Zk+ Zr's
Z(1,1) = (2*Zr + 2*Zk + Zm);
Z(1,2) = -(Zk+Zr);

if N > 2
    if N > 3
        for r = 2:N-1
            Z(r,r-1)= Z(1,2);
            Z(r,r)  = Zm-2*Z(1,2);
            Z(r,r+1)= Z(1,2);
        end
    end
    Z(N,N-1) = Z(1,2);
    Z(N,N) = Z(1,1);
end
Z = zpk(Z);
Y = inv(Z);
H = Y/s;
H = minreal(H);

%% ========== FREQUENCY RESPONSE SHAPES ================

% Do Nyquist
[Re,Im] = nyquist(H,omega);
 
%plot imaginary
figure("Name", "Imaginary Response", "NumberTitle", "off")
semilogx(omega,squeeze(Im(1,1,:))) %need to squeeze Im(1,1,:) into 2D
%%

%Magic kevin told me to do
H_w = Re + (1j)*Im; %Apparently this recombines the stuff
figure("Name", "Real Response", "NumberTitle", "off")
semilogx(omega,20*log10(squeeze(abs(H_w(1,1,:)))));
%findpeaks
[pks,index] = findpeaks(squeeze(abs(H_w(1,1,:))));
%pks are the values of the amplitudes, index will contain the indices
w_n = omega(index);

vv = imag(squeeze(H_w(1:length(index),1,:)).');
v = vv(index,:).';

[d1,d2,d3] = size(Im);
count = size(index);
u = zeros(d2,d1);
for X = 1:d1
    for Y = 1:count
        u(X,Y) = Im(X,1,index(Y));
    end
end

%%
figure(200)
bode(H(1,1,:));
%% Yeet 
% figure("Name", "Modes", "NumberTitle", "off")
% for mode=1:length(u)
%     subplot(2,2,mode)
%     plot([0; u(:,mode)], [0:4])
%     xline(0, "--")
%     ylim([0,4.5])
% end

% x1 = u1*exp(j*wn*t) or something like that

% Each mode crontributes a different amount to the system, so be able to
% alter where on the string you pluck.