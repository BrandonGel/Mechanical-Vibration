clc
clear all
close all
%%

m = 90;  %Table Mass
k = 10000; %Table Spring Constant
R = 10;
ma = 100; %Hanging Mass
ka = 10000; %Hanging Spring Constant
Ra = 10;

s = tf('s');

Z = [m*s + k/s + ka/s+R+Ra, -ka/s-Ra;
     -ka/s-Ra, ka/s-Ra+ma*s;];
Z = zpk(Z);

Y = inv(Z); %Response TF Matrix
H = Y/s;
H = minreal(H);

om = (1e-1:1e-1:1e3);
[Re,Im] = nyquist(H,om);
H_w = Re + (1j)*Im;

ind = 0;
wa = sqrt(ka/ma);

for indx = 1:2
   for indy = 1:2
       ind = ind + 1;
       subplot(2,2,ind)
       semilogx(om,20*log10(squeeze(abs(H_w(1,1,:)))),'linewidth',1);

   end
end
%%
F = [1; 0];
X = H*F;
t = (0:2*pi/(2.5*wa):100);
x = impulse(X,t);
figure(1)
plot(t,x)







