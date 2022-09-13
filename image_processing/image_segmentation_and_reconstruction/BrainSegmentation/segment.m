function segment()
% SSY186 Diagnostic Imaging 2016 
% Author: Your Name

close all; 
% define tissue labels
% see http://brainweb.bic.mni.mcgill.ca/brainweb/anatomic_normal.html
BCK=0; CSF=1;  GM=2; WM=3; FAT=4; MUSCLE=5; SKIN=6; SKULL=7; GLIA=8; CON=9;

% define filenames (using Brainweb names)
fnT1  = 't1_icbm_normal_1mm_pn0_rf0';   % T1
fnT2  = 't2_icbm_normal_1mm_pn0_rf0';   % T2
fnPD  = 'pd_icbm_normal_1mm_pn0_rf0';   % PD
fnGT  = 'phantom_1.0mm_normal_crisp';   % Ground Truth

% loadmnc from: http://www.mathworks.com/matlabcentral/fileexchange/32644-loadminc
% load data volumes
[T1,scaninfo] = loadminc([ fnT1 '.mnc']);
[T2,scaninfo] = loadminc([ fnT2 '.mnc']);
[PD,scaninfo] = loadminc([ fnPD '.mnc']);
[GT,scaninfo] = loadminc([ fnGT '.mnc']);

GT = uint8(GT);
GT(GT==GLIA) = GM;  % let's glia cells treat as gray matter
 
%T1 = PD;

% take a look at the volumes
T1= permute(T1, [1 2 3]); % look at a different view
%browseVolume(T1); 
%browseVolume(T2); 

%browseVolume(X);
% browseVolume(PD); browseVolume(GT);
%return  % remove it to run clustering

[M N K ] = size(T1);
Kc = round(K/2); % central index
 
ONESLICE=1;  % to make it easier to demonstrate we use ONE central (axial)slice only
if ONESLICE
    % Let's work with ONE slice only
    % display the central axial slice
    figure(10); T1 = T1(:,:,Kc);  imagesc(T1); colormap gray; title(['T1 slice#: ' num2str(Kc)]);
    T1n = T1 - min(T1(:)); T1n = T1n/max(T1n(:)); 
    %imwrite(T1n, 'T1.jpg','jpg','Quality', 100);
    
    figure(20); T2 = T2(:,:,Kc);  imagesc(T2); colormap gray; title(['T2 slice#: ' num2str(Kc)]);
    T2n = T2 - min(T2(:)); T2n = T2n/max(T2n(:)); 
    %imwrite(T2n, 'T2.jpg','jpg','Quality', 100);
    
    figure(30); PD = PD(:,:,Kc);  imagesc(PD); colormap gray; title(['PD slice#: ' num2str(Kc)]); 
    PDn = PD - min(PD(:)); PDn = PDn/max(PDn(:)); 
    %imwrite(PDn, 'PD.jpg','jpg','Quality', 100);
     
    figure(40); GT = GT(:,:,Kc);  imagesc(GT); title(['GT, slice #:' num2str(Kc)]);
    
end
 

sizeIMG = size(T1); % remember size of the volume
T1 = T1(:); T2 = T2(:); PD = PD(:);  GT = GT(:); % make data one dimensional == easier to cluster

% construct label masks
BRAIN_MASK = ((GT==CSF) | (GT==GM) | (GT==WM));
CSF_MASK   = GT==CSF;
GM_MASK    = GT==GM;
WM_MASK    = GT==WM;

% we are interested in brain tissues only so
% let's treat non-brain voxels as NaN
T1(~BRAIN_MASK) = NaN;   
T2(~BRAIN_MASK) = NaN;
GT(~BRAIN_MASK) = NaN;
 
% normalize to 0..1, why?  otherwise some modality 
% may dominate distance measurements
T1 = T1-min(T1(:)); T1 = T1/max(T1(:));
T2 = T2-min(T2(:)); T2 = T2/max(T2(:));
PD = PD-min(PD(:)); PD = PD/max(PD(:)); 

%% take a look how data look like, use scatterplots
figure; set(gca,'FontSize', 14); 
plot(T1(BRAIN_MASK), T2(BRAIN_MASK),'kx', 'LineWidth', 2); 
grid on; title('Brain data: data to segment');
xlabel('T1'); ylabel('T2');
 

% take a look how data look like, use scatterplots
% this time using true labels
figure; set(gca,'FontSize', 14); 
plot(T1(CSF_MASK), T2(CSF_MASK),'rx', 'LineWidth', 2); hold on;
plot(T1(GM_MASK),  T2(GM_MASK), 'go', 'LineWidth', 2); 
plot(T1(WM_MASK),  T2(WM_MASK), 'bs', 'LineWidth', 2); 
legend('CSF', 'GM', 'WM');
grid on; title('Brain data: true labels');
xlabel('T1'); ylabel('T2');

 
% put data into one vector/matrix
X = [T1 T2];  % for segmentation using T1 and T2
X = [T1 T2 PD];

% 
nrClasses = 3; % desired # of output segments
% do segmentation using k-means algorithm
% [0 1; 0.5 0.5; 1 0] [0.18 0.98; 0.47 0.42; 0.63 0.27]
[cidx, ctrs] = kmeans(X, nrClasses, 'Start', [0.18 0.98 1; 0.47 0.42 0.87; 0.63 0.27 0.73]);  % by clustering !
cidx = uint8(cidx);

% show k-means result as scatterplot
figure; set(gca,'FontSize', 14); 
plot(X(cidx==CSF,1),X(cidx==CSF,2),'r.', ...
     X(cidx==GM,1) ,X(cidx==GM,2) ,'g.', ...
     X(cidx==WM,1) ,X(cidx==WM,2) ,'b.');
xlabel('T1'); ylabel('T2'); legend('CSF','GM', 'WM');         
title('Segmented data'); grid on;

% calculate confusion matrix
cm = confusionmat(cidx(:), GT(:))  
di = diceIndex(cm);  % calculate Dice index
disp('Modality = (T1, T2), Protocol=ICBM, Phantom_name=normal, Slice_thickness=1mm, Noise=0%, INU=0%');
disp(['Dice index :   CSF  GM   WM : ' num2str(di(2:4))]); 

 
T1   = reshape(T1, sizeIMG); % make 2D image again
GT   = reshape(GT, sizeIMG);
cidx = reshape(cidx, sizeIMG);
errors =(GT~=cidx);
brainPixels = (GT>0);
nrErrors = sum(errors(:));
nrBrainPixels = sum(brainPixels(:));
ratio = nrErrors/nrBrainPixels;
disp(['Total # of errors: ' num2str(sum(errors(:)))  ' (' num2str(100*ratio) '%)']);

if ONESLICE
   
    % visualize  segmentation results in spatial domain
    % plot erroneously assigned pixels on the original images
    figure; set(gca,'FontSize', 12); 
    imagesc([cidx errors;
             cidx GT]); axis off;
    title(['UPPER-LEFT : errors on segmented image  \newline' ...
           'UPPER-RIGHT: error pixels \newline' ...
           'LOWER-LEFT : segmentation result \newline' ...
           'LOWER-RIGHT: Ground Truth image']);
    hold on;
    [x, y] = find(errors);
    plot(y, x, 'wo');
     
    figure; % plot errors on T1 image
    imagesc(T1); hold on; plot(y, x, 'ro'); colormap gray;
    title(['Wrongly classified pixels in red \newline' ...
           'Dice index : CSF  GM   WM : ' num2str(di(2:4)) '\newline' ...
           'Total # of errors: ' num2str(nrErrors)]); 
    
    %for i=1:length(x)  text(y(i),x(i), num2str(i),'Color', [1 1 0]); end
end
return
 
 
function browseVolume(V)
    figure;
    [M N K ] = size(V);
    for k=1:K
        slice = V(:,:,k); imagesc(slice); colormap gray; title(k);  drawnow;
    end
return
