clc
clear all
close all
%%
s = tf('s');

m1 = 4000;
m2 = m1;
m3 = m1;
m4 = m1;

k1 = 5000;
k2 = k1;
k3 = k1;
k4 = k1;

R = 10;
Z = [ k1/s+m1*s+k2/s+2*R, -k2/s-R, 0 , 0;
    -k2/s-R, k2/s+m2*s+k3/s+2*R, -k3/s-R, 0;
    0, -k3/s-R, k3/s+m3*s+k4/s+2*R, -k4/s-R;
    0, 0, -k4/s-R, k4/s+m4*s+R;];
Y = inv(Z);
H = Y/s;


w = [0.1: 0.001: 10];
figure(1);
bode(H,w);
[reals, imags] = nyquist(H,w);
H_w = reals + 1j*imags;

[~,mode_inds] = findpeaks(squeeze(abs(H_w(1,1,:))));
modes = w(mode_inds);

%%
v = imag(squeeze(H_w(1:4,1,:)).');
v = v(mode_inds,:).';

[d1, d2, d3] = size(imags);

u = zeros(d2,d1);
for indx = 1:d1
    for indy = 1:d2
        u(indx,indy) = imags(indx,1,mode_inds(indy));
    end
end

figure(7)
%plot(v+[-1e-9,-1e-9,-1e-9,-1e-9;-5e-10,-5e-10,-5e-10,-5e-10; 5e-10,5e-10,5e-10,5e-10;1e-9,1e-9,1e-9,1e-9],[1,4,7,10;1,4,7,10;1,4,7,10;1,4,7,10;],'LineStyle', 'none', 'Marker', 'o');
%plot(v+[-10,-10,-10,10;-5,-5,-5,-5; 5,5,5,5;10,10,10,10],[1,4,7,10;1,4,7,10;1,4,7,10;1,4,7,10;],'LineStyle', 'none', 'Marker', 'o');
%xlim([-2e-9,2e-9]);
%xlim([20,20]);
ylim([0,12]);

figure(8);
for mode = 1:length(u)
    %subplot(2,2,mode)
    plot([0; v(:,mode)],[0:4]);
    hold on
    xline(0,"--");
end 
%%
figure(9);
plot([0,0,0,0;v],[0,1,2,3,4]);

figure(10);
t = (0:0.01:10);
for ii = 1:1:length(t)
    plot(real([0,0,0,0;v.*exp(1j*modes*t(ii))]),[0,1,2,3,4]);
   %plot(real([0,0,0,0;v.*cos(1j*modes*t(ii))]),[0,1,2,3,4]);
    xlim([-0.1,0.1]);
    pause(0.001);
end

%% Building Problem (textbook 4.4.3)
m1 = 4000;          % kg
m2 = 4000;          % kg
m3 = 4000;          % kg
m4 = 4000;          % kg

k1 = 5000;          % N/m
k2 = 5000;          % N/m
k3 = 5000;          % N/m
k4 = 5000;          % N/m

R = 10;

Z = [k1/s + m1*s + k2/s + R, -k2/s - R, 0, 0;
    -k2/s - R, k2/s + m2*s + k3/s + R, -k3/s - R, 0;
    0, -k3/s - R, k3/s + m3*s + k4/s + R, -k4/s - R;
    0, 0, -k4/s - R, k4/s + m4*s + R];

Y = inv(Z);

H = Y/s;

w = [0.1:0.001:100];
% w = [1:0.001:10];
[reals, imags] = nyquist(H, w);
H_w = reals + (1j)*imags;
% figure()
% bode(H, w)

% Find natural/resonant frequencies (modes)
[peaks, mode_indxs] = findpeaks(squeeze(abs(H_w(1,1,:))));
modes = w(mode_indxs);

[d1, d2, d3] = size(imags);
u = zeros(d2, d1);  % rows = response, columns = masses
% get amplitude vectors
for indx_x=1:d1
    for indx_y=1:d2
        u(indx_x,indx_y) = imags(indx_x,1,mode_indxs(indx_y));
    end
end
