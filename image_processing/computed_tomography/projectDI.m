clear all; close all; clc;
load('data.mat');
%% PART A - Back Projection

% A2. Compute and plot the image back projection of the sinogram data
% Creating a back projection
%  backProjection(g,Theta) - implements equation (6.13)
%  NOTE: the angle is in degrees
theta = 30;
gBack = backProjection(g,theta);
imagesc(gBack);
colormap gray
title('Backprojection at 30 degrees');

% A3. Compute and plot the back projection summation image.
% The following is code for for construction of backprojection
% summation images. We use different N values: [6 12 30 60 180].
% The loops implement equation (6.14) in the textbook.
figure
% 6 back projections
N = 6;
b = 0;
for i=30:30:180
  b = b + backProjection(g,i);
end
gBack = pi/(2*N)*b;
subplot(2,2,1);
imagesc(gBack);
colormap gray
title('6 backprojections');

% 12 back projections
N = 12;
b = 0;
for i=15:15:180
  b = b + backProjection(g,i);
end
gBack = pi/(2*N)*b;
subplot(2,2,2);
imagesc(gBack);
colormap gray
title('12 backprojections');

% 30 back projections
N = 30;
b = 0;
for i=6:6:180
  b = b + backProjection(g,i);
end
gBack = pi/(2*N)*b;
subplot(2,2,3);
imagesc(gBack);
colormap gray
title('30 backprojections');

% 60 back projections
N = 60;
b = 0;
for i=3:3:180
  b = b + backProjection(g,i);
end
gBack = pi/(2*N)*b;
subplot(2,2,4);
imagesc(gBack);
colormap gray
title('60 backprojections');

figure
% 180 back projections
N = 180;
b = zeros(100);
for i=1:1:180
  b = b + backProjection(g,i);
end
gBack = pi/(2*N)*b;
imagesc(gBack);
colormap gray
title('180 backprojections');

%% PART B - Filtered Back Projection (FBP)
% Looking at equation (6.23) the filtering term is implemented by a ramp
% (1) filter that is first created.
% (2) Next G is aquired from g by taking its Fourier Transform
%      Note: fftshift is used to centre result so the filter can work
%            correctly
% (3) A pointwise multiplication between: filter and G is done
% (4) The inverse FFT is taken
% (5) Then a summation over all the angles from 0 to pi is done
% In all we've implemented equation (6.23)

%creating a ramp filter |p|:
filter = zeros(144,1);
for k=1:72
  filter(k) = 73 - k;
end
for k=1:72
  filter(k+72) = 0.1 + k;
end
ramp = filter/71;
ramp(ramp==0)=0.001; % it is customary for the origin of the filter to be nonzero

%creating a hamming filter (for B3)
ham = hamming(144);

figure
subplot(1,2,1);
plot(ramp);
axis([1 144 0 1]);
title('Ramp filter');
subplot(1,2,2);
plot(ham);
axis([1 144 0 1]);
title('Hamming filter');

% 1D fourier Transform
G = zeros(144,180);
for i=1:180
  G(:,i) = fftshift(fft(g(:,i)));
end

%  Filter: inverse 1-D Fourier Transform, and multiplication with filter
%  coefficients. Taking the absolute value of the ifft due to imaginary
%  components.
F = zeros(144,180);
for i=1:180
  F(:,i) = real(ifft(ifftshift(ramp.*G(:,i)))); % magnitude
  Fham(:,i) = real(ifft(ifftshift((ham.*ramp).*G(:,i)))); % magnitude
end

% Backprojection and summation

figure
% 6 back projections
N = 6;
f = 0;
fham = 0;
for i=30:30:180
  f = f + backProjection(F,i);
  fham = fham + backProjection(Fham,i);
end
fBack = pi/(2*N)*f;
fhamBack = pi/(2*N)*fham;
subplot(4,2,1);
imagesc(fBack);
title('Ramp filtering (N=6)');
subplot(4,2,2);
imagesc(fhamBack);
title('Hamming filtering (N=6)');
colormap gray

% 12 back projections
N = 12;
f = 0;
fham = 0;
for i=15:15:180
  f = f + backProjection(F,i);
  fham = fham + backProjection(Fham,i);
end
fBack = pi/(2*N)*f;
fhamBack = pi/(2*N)*fham;
subplot(4,2,3);
imagesc(fBack);
title('Ramp filtering (N=12)');
subplot(4,2,4);
imagesc(fhamBack);
title('Hamming filtering (N=12)');
colormap gray

% 30 back projections
N = 30;
f = 0;
fham = 0;
for i=6:6:180
  f = f + backProjection(F,i);
  fham = fham + backProjection(Fham,i);
end
fBack = pi/(2*N)*f;
fhamBack = pi/(2*N)*fham;
subplot(4,2,5);
imagesc(fBack);
title('Ramp filtering (N=30)');
subplot(4,2,6);
imagesc(fhamBack);
title('Hamming filtering (N=30)');
colormap gray

% 60 back projections
N = 60;
f = 0;
fham = 0;
for i=3:3:180
  f = f + backProjection(F,i);
  fham = fham + backProjection(Fham,i);
end
fBack = pi/(2*N)*f;
fhamBack = pi/(2*N)*fham;
subplot(4,2,7);
imagesc(fBack);
title('Ramp filtering (N=60)');
subplot(4,2,8);
imagesc(fhamBack);
title('Hamming filtering (N=60)');
colormap gray


% 1D fourier Transform
G = zeros(144,180);
for i=1:180
  G(:,i) = fftshift(fft(g(:,i)));
end

%  Filter: inverse 1-D Fourier Transform, and multiplication with filter
%  coefficients. Taking the absolute value of the ifft due to imaginary
%  components.
F = zeros(144,180);
Fham = zeros(144,180);
for i=1:180
  F(:,i) = real(ifft(ifftshift(ramp.*G(:,i)))); % magnitude
  Fham(:,i) = real(ifft(ifftshift((ham.*ramp).*G(:,i)))); % magnitude
end

% Backprojection and summation
f = 0;
fham = 0;
for theta = 1:180
  f = f + backProjection(F,theta);
  fham = fham + backProjection(Fham,theta);
end

fBack = pi/(2*N)*f;
fhamBack = pi/(2*N)*fham;
figure
subplot(1,2,1);
imagesc(fBack);
title('Filtered Backprojection (Ramp)');
subplot(1,2,2);
imagesc(fhamBack);
colormap gray
title('Filtered Backprojection (Hamming)');

%% PART C - Convolution Back Projection (CBP)
% Here we are implementing equation (6.24) with a windowed filter from
% (6.27).
% (1) Take out of the zero value in the ramp filter (NOT USED!)
% (2) Take fourier transform of the window (Hamming window used)
% (3) perform the convolution, with the 'same' argument which will
%     extract the middle part of the convolution and drop the transients
% (4) Perform the summation of the backprojection.

% Convolution
figure
c_tilde = fftshift(ifft(ifftshift(ham.*ramp)));
plot(real(c_tilde));
title('Inverse Fourier Transform of the Hamming.*Ramp filter');

f = zeros(144,180);
for i=1:180
  f(:,i) = real(conv(c_tilde,g(:,i),'same'));
end

% Backprojection and summation
b = 0;
for i=1:180
  b = b + backProjection(f,i);
end

figure
fBack = pi/(2*N)*b;
subplot(1,2,2);
imagesc(fBack);
colormap gray
title('Convolution Backprojection (Hamming)');
subplot(1,2,1);
imagesc(fhamBack);
colormap gray
title('Filtered Backprojection (Hamming)');

%% E. Real CT-image
% Perform
%  (1) Back projection
%  (2) Filtered Bakcprojection
%  (3) Convolutional Backprojection
% on a real CT-image.
load('data2.mat');
figure
% {A} Back Projection
subplot(1,3,1);
N = 180;
b = 0;
for i=1:1:180
  b = b + backProjection2(g2,i);
end
gBack = pi/(2*N)*b;
imagesc(gBack);
colormap gray
title('Backprojection');

%creating the filter
filter = zeros(1095,1);
for k=1:547
  filter(k) = 548 - k;
end
for k=1:547
  filter(k+548) = 0.1 + k;
end
ramp = filter/548;
ramp(ramp==0)=0.001; % it is customary for the origin of the filter to be nonzero
ham = hamming(1095);
filter = ham.*ramp;
filter = (filter - min(filter)) / ( max(filter) - min(filter) ); % nromalizing between 0 and 1

% {B} Filtered Back Projection
subplot(1,3,2);
tic
% 1D fourier Transform
Gham = zeros(1095,180);
for i=1:180
  Gham(:,i) = fftshift(fft(g2(:,i)));
end
% Filter
Fham = zeros(1095,180);
for i=1:180
  Fham(:,i) = real(ifft(ifftshift(filter.*Gham(:,i)))); % magnitude
end
% Backprojection and summation
fham = 0;
for theta = 1:180
  fham = fham + backProjection2(Fham,theta);
end
fhamBack = pi/(2*N)*fham;
toc
imagesc(fhamBack);
colormap gray
title('Filtered Backprojection (Hamming)');

% {C} Convolution Back Projection
subplot(1,3,3);
tic
c_tilde = fftshift(ifft(ifftshift(filter)));
plot(real(c_tilde));
title('Inverse Fourier Transform of the Hamming.*Ramp filter');

f = zeros(1095,180);
for i=1:180
  f(:,i) = real(conv(c_tilde,g2(:,i),'same'));
end
% Backprojection and summation
b = 0;
for i=1:180
  b = b + backProjection2(f,i);
end
fBack = pi/(2*N)*b;
toc
imagesc(fBack);
colormap gray
title('Convolution Backprojection (Hamming)');

figure
plot(abs(filter));
title('Hammming.*Ramp filter, in Fourier space');