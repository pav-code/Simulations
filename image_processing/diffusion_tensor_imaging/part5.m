% Part 5: Analysis of MRI Data - DTI
clear all; close all; clc;
%1) Find FA and PDD of diffusion tensor (D)
D = 10^-4*[17 0 0; 0 7 0; 0 0 4];
[eVec, eVal] = eig(D); % calculate eigenvalues and eigenvectors

% Fractional Anisotropy
FA = sqrt(3/2)*sqrt(sum(sum((eVal-mean(mean(eVal)))^2))/(sum(sum(eVal^2))));

% Principal Direction of Diffusion
PDD = eVec(:,1);
angErr = zeros(100,1);
PDD_hat = zeros(3,100);
FA_hat = zeros(100,1);

g = random_unit_vector(30); % MATLAB function 
for l = 1:100
%2) Simulate the diffusion signal arising from D

b = 1500; % diffusion weighting
s0 = 1000; % signal without diffusion sensitizing gradient

% Gaussian noise
SNR = 15; 
var = s0/SNR;
n1 = var*randn(30,1);
n2 = var*randn(30,1);

% Diffusion signal (S) measured in direction g
s = zeros(30,1);
for i=1:30
    s(i) = s0*exp(-b*g(:,i)'*D*g(:,i));
end

% Noisy signal (Sm)
Sm = sqrt((s+n1).^2+n2.^2);

% LS estimation of D_hat
y = -log(Sm/s0)/b; % measurements
G = zeros(30,6);
for i=1:30
    G(i,:) = [g(1,i)^2 g(2,i)^2 g(3,i)^2 ...
        2*g(1,i)*g(2,i) 2*g(1,i)*g(3,i) 2*g(2,i)*g(3,i)];
end
t = G\y;
D_hat = [t(1) t(4) t(5); 0 t(2) t(6); 0 0 t(3)];
[eVec, eVal] = eig(D_hat); % calculate eigenvalues and eigenvectors

% Fractional Anisotropy of the estimated D
FA_hat(l) = sqrt(3/2)*sqrt(sum(sum((eVal-mean(mean(eVal)))^2))/(sum(sum(eVal^2))));

% Principal Direction of Diffusion of the estimated D
PDD_hat(:,l) = eVec(:,1);

% Angular error
angErr(l) = rad2deg(atan2(norm(cross(PDD,PDD_hat(:,l))),dot(PDD,PDD_hat(:,l))));end
