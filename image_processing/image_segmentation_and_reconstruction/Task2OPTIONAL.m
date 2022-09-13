close all; clear all; clc;

N = 100;          % The discretization points.
theta = 0:1:179;  % No. of used angles.
p = 75;           % No. of parallel rays.
eta = 0.05;       % Relative noise level.
k = 100;           % No. of iterations.

phantom_rgb = imread('phantom.png');
phantom = rgb2gray(phantom_rgb)/255;

[A,b_ex,x_ex] = paralleltomoOPTIONAL(N,theta,p);
e = 0; % 0 error
b = b_ex + e;

% Show the exact solution.
figure
imagesc(reshape(x_ex,N,N)), colormap gray,
axis image off
c = caxis;
title('Exact phantom')

fprintf(1,'\n\n');
fprintf(1,'Perform k = %2.0f iterations with Kaczmarz''s method.',k);
fprintf(1,'\nThis takes a moment ...');
% Perform the kaczmarz iterations.
Xkacz = kaczmarz(A,b,k);

% Show the kaczmarz solution.
figure
imagesc(reshape(Xkacz,N,N)), colormap gray,
axis image off
title('Kaczmarz reconstruction')