close all;
clear;
clc;

%% Specifications
wp1  =   0.3*pi;
wp2  =   0.6*pi;
wa1  =   0.1*pi;
wa2 =    0.8*pi;

n1  =   12;
n2  =   12;
mu  =   0.1;
epsilon_t   =   4.42e-4;
w   =   1;
M   =   100*100;
m   =   sqrt(M);
plot_samples    = 128;

%% Define parameters

n       =   (n2+1)*(2*n1+1);
f       =   [1;zeros(n,1)];

[wx,wu] =   meshgrid(-pi:2*pi/(m-1):pi,-pi:2*pi/(m-1):pi);

alpha1 = 70;
alpha2 = 20;
theta = 5;
T = 0.02*pi;
B = 0.9*pi;
Ba = 0.96*pi;

m1 = tand(alpha1-theta);
m2 = tand(alpha1+theta);
c1 = T/cosd(alpha1-theta);
c2 = T/cosd(alpha1+theta);
ca1 = (T+0.07*pi)/cosd(alpha1-theta);
ca2 = (T+0.07*pi)/cosd(alpha1+theta);

md1 = tand(alpha2-theta);
md2 = tand(alpha2+theta);
cd1 = T/cosd(alpha2-theta);
cd2 = T/cosd(alpha2+theta);
cda1 = (T+0.07*pi)/cosd(alpha2-theta);
cda2 = (T+0.07*pi)/cosd(alpha2+theta);

Ad = ((wu>=m1*wx-c1 & wu<=m2*wx+c2)|(wu<=m1*wx+c1 & wu>=m2*wx-c2))&(wu.^2+wx.^2<=B^2)...
    |((wu>=md1*wx-cd1 & wu<=md2*wx+cd2)|(wu<=md1*wx+cd1 & wu>=md2*wx-cd2))&(wu.^2+wx.^2<=B^2);
Ad_s = ((wu>=m1*wx-ca1 & wu<=m2*wx+ca2)|(wu<=m1*wx+ca1 & wu>=m2*wx-ca2))&(wu.^2+wx.^2<=Ba^2)...
    |((wu>=md1*wx-cda1 & wu<=md2*wx+cda2)|(wu<=md1*wx+cda1 & wu>=md2*wx-cda2))&(wu.^2+wx.^2<=Ba^2);
W       =  (1-(Ad_s))*w+Ad;

% wp  =   0.5*pi;
% wa  =   0.7*pi;
% W       =   ((w1.^2 + w2.^2) <= wp^2 ) + (((w1.^2 + w2.^2) >= wa^2 )*w);
% Ad      =   (w1.^2 + w2.^2) <= wp^2;

contour(wx,wu,Ad)
hold on;
contour(wx,wu,W)
grid on
%% 
dropout =   sum(sum(W==0));


count_i =1;
M_act       =  M-dropout;
A1_temp     = zeros(M_act,n);
b1_temp     = zeros(M_act,1);
for p = 1:M
    w1_p    = wx(p); 
    w2_p    = wu(p);
   
    n2_array = 0:n2;
    cos_w1w2 = [];
    W_k = (W(p));
    if W_k
        for m1= -n1:n1
            theta = w2_p*n2_array + m1*w1_p;
            cos_theta = cos(theta);
            cos_w1w2 = W(p)*[cos_theta cos_w1w2];
        end

        A1_temp(count_i,:)  = W(p)*cos_w1w2;
        b1_temp(count_i)    = W(p)*Ad(p);
        count_i = count_i +1; 
    end
end

A1          =   [ones(M_act,1) A1_temp;
                 ones(M_act,1) -A1_temp];
b1          =   [b1_temp; -b1_temp];


 
cvx_begin
    variable x(n+1)
    minimize (f'*x)
    subject to
        A1 * x  >= b1
cvx_end

h = x(2:end);
H = flipud(reshape(h,n2+1,2*n1+1));
H_hat= [H(1:n2,:)/2;H(n2+1,:);flipud(fliplr(H(1:n2,:)/2))];
save('H_hat.mat','H_hat')

figure;
surf(H_hat);

e_plot      = ones(plot_samples,1);
[Amp,fx,fy] = freqz2(H_hat,[plot_samples plot_samples]);
gain        = 20*log10(abs(Amp));
figure;
mesh(fx,fy,gain,gain.*(gain>-100))
zlim([-100 5]);
