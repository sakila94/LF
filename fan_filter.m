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
% W       =  ((w1 >= wp1 & w1 <= wp2 & w2 >= wp1 & w2 <= wp2)|(w1 <= -wp1 & w1 >= -wp2 & w2 <= -wp1 & w2 >= -wp2) +...
%              ((w1 <= wa1 | w1 >= wa2 | w2 <= wa1 | w2 >= wa2)&(w1 >= -wa1 | w1 <= -wa2 | w2 >= -wa1 | w2 <= -wa2))*w);
% 
% Ad      =   (w1 >= wp1 & w1 <= wp2 & w2 >= wp1 & w2 <= wp2)|(w1 <= -wp1 & w1 >= -wp2 & w2 <= -wp1 & w2 >= -wp2);

% mp1 = 1;
% mp2 = 0.5;
% ma1 = 0.3;
% ma2 = 1.5;
% wp  =   0.8*pi;
% wa = 0.9*pi;
% c= 0.2*pi;
alpha = 30;
theta = 10;
T = 0.08*pi;
B = 0.9*pi;
Ba = pi;

m1 = tand(alpha-theta);
m2 = tand(alpha+theta);
c1 = T/cosd(alpha-theta);
c2 = T/cosd(alpha+theta);
ca1 = (T+0.1*pi)/cosd(alpha-theta);
ca2 = (T+0.1*pi)/cosd(alpha+theta);

Ad = ((wu>=m1*wx-c1 & wu<=m2*wx+c2)|(wu<=m1*wx+c1 & wu>=m2*wx-c2))&(wu.^2+wx.^2<=B^2);
W       =  (1-( ((wu>=m1*wx-ca1 & wu<=m2*wx+ca2)|(wu<=m1*wx+ca1 & wu>=m2*wx-ca2))&(wu.^2+wx.^2<=Ba^2)))*w+Ad;

% wp  =   0.5*pi;
% wa  =   0.7*pi;
% W       =   ((w1.^2 + w2.^2) <= wp^2 ) + (((w1.^2 + w2.^2) >= wa^2 )*w);
% Ad      =   (w1.^2 + w2.^2) <= wp^2;

surf(wx,wu,Ad)
hold on;
surf(wx,wu,W)
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
