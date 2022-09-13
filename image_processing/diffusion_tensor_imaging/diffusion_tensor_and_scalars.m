clc
clear all
close all

% SET THE NUMBER OF DIRECTIONS USED IN THE DTI SCAN HERE!!!!
number_of_directions = 25; 

% This data originates from a breast DTI scan. There were 40 axial slices.
% The data were stored in the following DICOM files:
%
% DICOM file                                  Description
% ----------------------------------------------------------------------------
% IM-0009-0001.dcm .... IM-0009-0040.dcm      b=0 image for each slice
% IM-0009-0041.dcm .... IM-0009-0080.dcm      direction 1 image for each slice
% IM-0009-0081.dcm .... IM-0009-0120.dcm      direction 2 image for each slice
% IM-0009-0121.dcm .... IM-0009-0160.dcm      direction 3 image for each slice
%        .                     .              .
%        .                     .              .
%        .                     .              .
%        .                     .              .
%
% In this exercise you only have the data for slice 1 and slice 19.

number_of_axial_slices = 40; % number of axial slices acquired
which_slice = 19; % let's only consider slice 19

oldfolder = cd('breastDTIdata');

% Obtain the unit vector for each direction. We do this by examining the
% DICOM header for each of the files (except b=0) for slice 1. The
% information is stored in vendor-specific private fields.

dirs = zeros(number_of_directions, 3);

for i = 1:number_of_directions
   n = i * number_of_axial_slices + 1; % 41, 81, etc
   file_name = ['IM-0009-',num2str(n, '%0.4d'),'.dcm'];
   disp(file_name)
   info = dicominfo(file_name);
   dirs(i, 1) = info.Private_0019_10bb;
   dirs(i, 2) = info.Private_0019_10bc;
   dirs(i, 3) = info.Private_0019_10bd;
end


% INSERT YOUR CODE HERE TO SHOW THAT EACH DIRECTION IS A UNIT VECTOR
mag = zeros(25,1);
for i=1:number_of_directions
    mag(i) = sum(dirs(i,:).^2); % the magnitude of each gradient is approx. 1
end

% Determinbe the b value used. Again tbis is stored in a vendor-specific
% private field
b = info.Private_0043_1039(1);

% Create a matrix B such that each row contains the 6 bij for a given
% direction. See Park equations 8 and 9 and slide 5 of Westin.

B = zeros(number_of_directions, 6);
for i = 1:number_of_directions
    g = dirs(i,:);
    gg = g' * g;
    B(i,1) = gg(1,1);
    B(i,2) = gg(2,2);
    B(i,3) = gg(3,3);
    B(i,4) = 2 * gg(1,2);
    B(i,5) = 2 * gg(1,3);
    B(i,6) = 2 * gg(2,3);
end

B = B * b;

% Read in the baseline image (b=0) for slice 19
file_name = ['IM-0009-',num2str(which_slice, '%0.4d'),'.dcm'];
disp(file_name)
base_img = double(dicomread(file_name));

% INSERT YOUR CODE HERE TO DISPLAY THE b=0 IMAGE FOR SLICE 19
imagesc(base_img);
colormap gray
title('b=0 image for the slice 19');

% S array contains all of the data (except for S0)
S_array = zeros([number_of_directions size(base_img)]);
for i = 1:number_of_directions
    n = i * number_of_axial_slices + 19;
    file_name = ['IM-0009-',num2str(n, '%0.4d'),'.dcm'];
    disp(file_name);
    S_array(i,:,:) = dicomread(file_name);
end


% Now compute the diffusion tensor for each voxel (see Le Bihan JMRI
% 2001, p535) and some scalar measures for this tensor.

mean_diffusivity = zeros(size(base_img)); % mean diffusivity
FA  = zeros(size(base_img)); % fractional anisotropy
VR  = zeros(size(base_img)); % volume ratio


for x = 1:size(base_img, 1)
    for y = 1:size(base_img, 2)
        
        S = S_array(:,x,y);
        S0 = base_img(x,y);

        % compute the diffusion tensor for this voxel
        if S0
            
            S = S ./ S0;

            D = -B \ log(S); %same as -inv(B) * log(S)

            DT = zeros(3);
            DT(1,1) = D(1);
            DT(2,2) = D(2);
            DT(3,3) = D(3);
            DT(2,1) = D(4);
            DT(1,2) = D(4);
            DT(1,3) = D(5);
            DT(3,1) = D(5);
            DT(2,3) = D(6);
            DT(3,2) = D(6);
            
            % compute some scalar parameters for the voxel
            
            if ~sum(isnan(DT(:))) %avoid cases where the tensor contains nans
              
                [evect, evals] = eig(DT); %find eigenvalues/vectors
                
                if (evals(1,1) > 0 && evals(2,2) > 0 && evals(3,3) > 0) 
                    mean_diffusivity(x,y) = (evals(1,1) + evals(2,2) + evals(3,3))/3;
                    num = sum(eye(3).*(evals - mean_diffusivity(x,y)).^2);
                    den = sum(evals.^2);
                    
                    % INSERT YOUR CODE HERE TO COMPUTE: 
                    % (1) fractional anisotropy (FA), and
                    % (2) volume ratio (VR).
                    FA(x,y) = sqrt(3/2)*sqrt(num/den);  
                    VR(x,y) = (evals(1,1)*evals(2,2)*evals(3,3))/mean_diffusivity(x,y).^3;
                end          
            end
        end
    end
end
figure
subplot(1,3,1);
imshow(FA); title('Fractional Anisotropy (FA)');
subplot(1,3,2);
imshow(VR); title('Volume Ratio (VR)');
subplot(1,3,3);
imshow(mean_diffusivity, []); title('Mean Diffusivity (MD)');
cd(oldfolder)
