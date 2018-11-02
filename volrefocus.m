% function volrefocus

% VOLREFOCUS implements the 4-D fan filters for volumetric refocusing of
% LFs.
% The fan filters hxu, hxuth and light field LF are loaded from directory
% and converted to frequency domain using FFT. Filtering is done in 4-D
% frequncy domain to produce outputs lfvf with conventional depth filter
% and lfvfh with hard thresholded depth filter (sparse filter).
% Output images lfvf and lfvfh are saved in .mat format.
%
% Inputs:
%       lfnm - Light field name
%       hxu - Fan filter impulse response. Loading from saved directory.
%       hxuht - Hard thresholded fan filter impulse response. Loading from
%               saved directory
% 
% Outputs:
%       lfvf - Refocused image from conventional depth filter
%       lfvfh - Refocused image from hard thresholded depth
%               filter (sparse filter)
%
% Author - Chamira Edussooriya
% Edited by - Sanduni U. Premaratne
% Date - Apr 07, 2018
% Last modified - May 24, 2018

% Edited by - Sakila Jayaweera
% Date - October 25, 2018

clc;
close all;
clear all;

path = 'D:\SAKILA\Research\Multiple_reforcasing_LF\Light';
lfnm = [path,'\Books'];         % name of the LF
% lfnm = 'Sophie_&_Vincent_1';
load([lfnm,'.mat']);      % load the LF

LF = LF(3:13,3:13,:,:,1:3); % #ok<NODEF>
LF = im2uint8(LF);

% [Multi, MultiSp] = multi_fan_filter(30,-50,5);
[Multi, MultiSp] = fanFilterSparse(50,10);
% load('H_hat50.mat');    % load coefficients of conventional fan filter
% hxu = H_hat;
% 
% load('H_hat50n.mat');      % load coefficients of hardthresholed fan filter
% hxuht = H_hat;

hxu = Multi;
hxuht = MultiSp;
hyv = hxu;              % the same filter used for both Hxu and Hyv
hyvht = hxuht;          % the same filter used for both Hxuht and Hyvht


[Ny,Nx,Nv,Nu,Nc] = size(LF);
if Nc == 4              % remove the weight channel
    Nc = 3;
end
[Mx,Mu] = size(hxu);
[My,Mv] = size(hyv);

Fx = Mx + Nx - 1;           % length of the fftn for x dimension 
Fy = My + Ny - 1;           % length of the fftn for v dimension
Fu = Mu + Nu - 1;           % length of the fftn for u dimension 
Fv = Mv + Nv - 1;           % length of the fftn for v dimension

if max(LF(:))<256           % tesnor to store volumetric focussed LFs
    lfvf = zeros(Nv,Nu,Nc,'uint8');  
    lfvfh = zeros(Nv,Nu,Nc,'uint8');  
else
    lfvf = zeros(Nv,Nu,Nc,'uint16');
    lfvfh = zeros(Nv,Nu,Nc,'uint16');
end

disp('Implement Conventional 4-D Fan Filter'); 

Hxu = fft2(hxu,Fx,Fu);
Hyv = fft2(hyv,Fy,Fv);
Hxu = permute(repmat(Hxu,[1,1,Ny,Nv]),[3,1,4,2]);  
Hyv = repmat(Hyv,[1,1,Nu]);  
    
for c = 1:Nc  
    disp(['Vol Refocusing, c = ',num2str(c)]);      
        
    lfin = squeeze(LF(:,:,:,:,c));
    LFIN = fft(fft(lfin,Fx,2),Fu,4);
    clear('lfin');
    LFOUT = LFIN .* Hxu;
    clear('LFIN');
    lfout = ifft(ifft(LFOUT,Fx,2),Fu,4);
    lfout = squeeze(lfout(:,(Fx+1)/2,:,(Mu+1)/2:end-(Mu-1)/2));
    clear('LFOUT');
    
    LFIN = fft(fft(lfout,Fy,1),Fv,2);
    clear('lfout');
    LFOUT = LFIN .* Hyv;
    clear('LFIN');
    lfout = ifft(ifft(LFOUT,Fy,1),Fv,2);
    clear('LFOUT'); 
    lfvf(:,:,c) = real(squeeze(lfout((Fy+1)/2,(Mv+1)/2:end-(Mv-1)/2,:)));
    clear('lfout');
end
save([lfnm,'_VF.mat'],'lfvf');
imshow(lfvf);
clear('lfvf');

disp('Implement Hardthresholded 4-D Fan Filter');  

Hxu = fft2(hxuht,Fx,Fu);
Hyv = fft2(hyvht,Fy,Fv);
Hxu = permute(repmat(Hxu,[1,1,Ny,Nv]),[3,1,4,2]);  
Hyv = repmat(Hyv,[1,1,Nu]);  
    
for c = 1:Nc  
    disp(['Vol Refocusing, c = ',num2str(c)]);      
        
    lfin = squeeze(LF(:,:,:,:,c));
    LFIN = fft(fft(lfin,Fx,2),Fu,4);
    clear('lfin');
    LFOUT = LFIN .* Hxu;
    clear('LFIN');
    lfout = ifft(ifft(LFOUT,Fx,2),Fu,4);
    lfout = squeeze(lfout(:,(Fx+1)/2,:,(Mu+1)/2:end-(Mu-1)/2));
    clear('LFOUT');
    
    LFIN = fft(fft(lfout,Fy,1),Fv,2);
    clear('lfout');
    LFOUT = LFIN .* Hyv;
    clear('LFIN');
    lfout = ifft(ifft(LFOUT,Fy,1),Fv,2);
    clear('LFOUT'); 
    lfvfh(:,:,c) = real(squeeze(lfout((Fy+1)/2,(Mv+1)/2:end-(Mv-1)/2,:)));
    clear('lfout');
end
save([lfnm,'_VFHT.mat'],'lfvfh');
figure, imshow(lfvfh);
clear('lfvf');

