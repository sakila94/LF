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
e_n     =   ones(n,1);
f       =   [1; mu*e_n; zeros(n,1)];

[wx,wu] =   meshgrid(-pi:2*pi/(m-1):pi,-pi:2*pi/(m-1):pi);

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

A1          =   [ones(M_act,1) zeros(M_act,n) A1_temp;
                 ones(M_act,1) zeros(M_act,n) -A1_temp];
b1          =   [b1_temp; -b1_temp];

A2          =   [zeros(n,1) eye(n) eye(n);
                zeros(n,1) eye(n) -eye(n)];
            
A           =   [A1; A2];
b           =   [b1; zeros(2*n,1)];

 
cvx_begin
    variable x(2*n+1)
    minimize (f'*x)
    subject to
        A * x  >= b
cvx_end

h = x(n+2:2*n+1);
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
%%

H_ht        = zeros(n2+1,2*n1+1);
H_ht(1:n2,:) = H(1:n2,:)/2.* ( abs(H(1:n2,:)/2) >= epsilon_t);
H_ht(n2+1,:) = H(n2+1,:).* ( abs( H(n2+1,:)) >= epsilon_t);
% H_ht(1,:)   = [h22 h23_T].* ( abs([h22 h23_T]) >= epsilon_t);
% H_ht(:,1)   = [h22; h32].* ( abs([h22; h32]) >= epsilon_t);
% H_ht(2:end,2:end)   = H33 .* ( abs(H33) >= epsilon_t);
H_ht = flipud(H_ht);
% counting zeros
count       = sum(sum(H_ht(2:end,2:end) == 0));
count       = (sum(sum(H_ht == 0)) - count) *2 + count*4
%% Calculate required parameters

h_ht        = H_ht(:);
tau         = find(h_ht==0);
tau_bar     = find(h_ht~=0);
n_sp        = length(tau_bar)  %length(h_ht)  - length(tau)

for h_i = flipud(tau)
    A1_temp(:,h_i) = [];
end

A_new   =   [ones(M_act,1) A1_temp;
             ones(M_act,1) -A1_temp];
b_new   =   b1; 
f_new   =   [1; zeros(n_sp,1)];

cvx_begin
    variable x_new(n_sp+1)
    minimize (f_new'*x_new)
    subject to
        A_new * x_new  >= b_new
cvx_end

%% Reconstructing the impulse response

h_new   = zeros(n,1);
h_new(tau_bar) = x_new(2:end);

H_new = flipud(reshape(h_new,n2+1,2*n1+1));
H_hatSp= [H_new(1:n2,:)/2;H_new(n2+1,:);flipud(fliplr(H_new(1:n2,:)/2))];
save('H_hatSparse.mat','H_hatSp')
         
%% Visualization of results

figure;
surf(H_hatSp);
         
[Amp,fx,fy] = freqz2(H_hatSp,[plot_samples plot_samples]);
gain        = 20*log10(abs(Amp));

figure;
mesh(fx,fy,gain,gain.*(gain>-80))
zlim([-80 5]);
