clc
close all;
clear all;
width=550;
height=400;

%% Building Problem (textbook 4.4.3)
m1 = 4000;          % kg
m2 = 4000;          % kg
m3 = 4000;          % kg
m4 = 4000;          % kg

k1 = 5000;          % N/m
k2 = 5000;          % N/m
k3 = 5000;          % N/m
k4 = 5000;          % N/m

s = tf('s');

R = 100;
Z = [ (k1+k2)/s+m1*s+2*R, -k2/s-R, 0 , 0;
    -k2/s-R, (k2+k3)/s+m2*s+2*R, -k3/s-R, 0;
    0, -k3/s-R, (k3+k4)/s+m3*s+2*R, -k4/s-R;
    0, 0, -k4/s-R, k4/s+m4*s+R;];

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


v = imag(squeeze(H_w(1:4,1,:)).');
v = v(mode_indxs,:).';

%%
tt = (0:0.01:10);
for ii = 1:length(tt)
    figure(8);
    for mode = 1:length(v)-3
        %subplot(2,2,mode)
        plot([0; v(:,mode)*cos(tt(ii))],[0:4]);
        xline(0,"--");
        xlim([-0.08,0.08]);
    end 
end
%%
x0 = [0.001; 0.01; 0.02; 0.025;];
X0 = [-0.0107160745085793;-0.0201380441946668;-0.0271311774098532;-0.0308521127434055];
x00 = [x0;0];
X0 = x0/s;
k = [k1,k2,k3,k4];
% for ii = 1:length(x0)
%     F(ii,1) = -k(ii)*(x00(ii+1)-x00(ii))/s;
% en

f1= -k1*x0(1);
f2 = -k2*(x0(2)-x0(1));
f3 = -k3*(x0(3)-x0(2));
f4 = -k4*(x0(4)-x0(3));
F(1,1) = f1-f2;
F(2,1) = f2-f3;
F(3,1) = f3-f4;
F(4,1) = f4;
F = F/s;
X = H*F +X0;
t = (0:0.1:1000);
x = impulse(X,t);

figure(2);
plot(t,x);

%%
for ii = 1:length(t)
    figure(3);
    plot([0,x(ii,:)],[0,1,2,3,4]);
    xlim([-0.08,0.08]);
    pause(0.01);
end
%% %%%%%%%%%%%%Rotating Spring Rod%%%%%%%%%%%%
k1 = 1000;
k2 = 1000;
m = 10;
R = 1;
s = tf('s');

Z = [ 1/3*m*s+k1/s+R, 1/6*m*s;
      1/6*m*s, 1/3*m*s + k2/s+R;];
Y = inv(Z);
H = Y/s;

w = [0.01:0.001:100];
[real,imags] = nyquist(H,w);
H_w = real+1j*imags;

[d1, d2, d3] = size(H_w);

figure("Name", "H_w(1,1) Real", "NumberTitle", "off")
semilogx(w,20*log10(abs(squeeze(H_w(1,1,:)))));
set(gcf,'position',[10,50,width,height]);

figure("Name", "H_w Real", "NumberTitle", "off")
ind = 1;
for x = 1:d1
    for y = 1:d2
        subplot(2,2,ind);
        semilogx(w,20*log10(abs(squeeze(H_w(x,y,:)))));
        ind = ind + 1;
    end
end
set(gcf,'position',[600,50,width,height]);

% Find natural/resonant frequencies (modes)
[peaks, mode_indxs] = findpeaks(squeeze(abs(H_w(1,1,:))));
modes = w(mode_indxs);

v = imag(squeeze(H_w(1:2,1,:)).');
v = v(mode_indxs,:).';

figure("Name", "Modes", "NumberTitle", "off")
for mode = 1:length(v);
   subplot(2,1,mode);
   plot([1:2],[v(:,mode)]);
   ylim([-0.05,0.05]);
end
set(gcf,'position',[10,540,width,height]);

%% Model Animation
frames = 1000;
time = linspace(0,10,frames);

figure("Name", "modes", "NumberTitle", "off");
for indx = 1:frames
    plot([1:2],[v.*cos(modes*time(indx))]);
    yline(0,"--");
    xlim([1,2]);
    legend(["Mode 1", "Frame: " + indx + "/" + frames])
    ylim([-.1, .1]);
    pause(1e-5);
end
    












