function BW = segmentDAPIimage(IM,ROI,approxPixelWidthOfCellNuclei,manual_input_yes_no,answers)
% segmentDAPIimage - returns binary mask of cell nuclei from input DAPI image (of type double, uint8 or uint6)
%
% INPUTS:
%
% IM: input image, DAPI channel 2D image - array of type double, uint8 or
% uint6. 
%
% ROI: Mask with ones everywhere where you want to segment nuclei (give emty image - all zeros - if you want nuclei segmented everywhere).
% approxPixelWidthOfCellNuclei: The average diameter measured in # of
% pixels of the DAPI stained nuclei you want to segment (approx 25 pixels for our test images)
%
% manual_input_yes_no: if set to one 1 you will be promtet 4 time during
% the run of the function and have to pick between images that show different
% theshold levels. (Eg. you will be shown a fig. and asked "What thresholded image includes all nuclei and not more? Note: Nuclei should 
% not have holes in them. (Enter image number from fig. here and press
% enter) :").
% if manual_input_yes_no is set to 0 it will run without promting using values optimized for our
% test images. 
% If you want to change these values in order to use this function for a large batch of your own images 
% without getting promtet, then search in this text for 'promt' and locate and change manually where: 
% opt_index = 3, pess_index = 9, regT_index = 7, altPess_index = 11).
% and set them to the values that you find fits your sets of images best, and theb run with manual_input_yes_no=0.
%
% answers: vector [opt_index pess_index regT_index altPess_index] with 4 integer index numbers in the range [1-20]. Each
% correspond to the appropriate manual command window replies given during
% a run with 'manual_input_yes_no=1'. Use this if you want to differentiate
% between different groups of images in a batch run. Note if set to [0 0 0 0]
% the default [3 9 7 10] will be used.
%
% OUTPUT:
% Binary image with 1's where there is nuclei and where 'touching' nuclei
% have been seperated by watershed
%
% NOTE: at the end of this function a plot of the final segmentation result
% is generated. Comment this out if you want to supress this during a batch
% run.

disp('.')
disp('Starting...')

opt_index     = answers(1);
pess_index    = answers(2);
regT_index    = answers(3);
altPess_index = answers(4);


IM = double(IM);% cast to double

approxPixelWidthOfCellNuclei = round(1.6*approxPixelWidthOfCellNuclei);

if sum(ROI(:))==0 % if this mask is empty make it an all ones image
    ROI = ones(size(IM));
end

% gauss and median filt DAPI channel image (to reduce noise on a scale much smaller than a nucleus)
IMg = medfilt2(imgaussfilt(IM,[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(0.75/50)),floor([approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(3/50)));% Find out if background levels are uneven
IMg(IMg<0) = 0;% gauss filt make some values go negative get rid of those


IMgvec = IMg(:);
IMgvec(ROI(:)==0)=NaN;% set pixel values outside epithlium to NaN

% normalize image with respect to max inside ROI
IMg = IMg./nanmax(IMgvec(:));

IMgvec = IMg(:);
IMgvec(ROI(:)==0) = NaN;% set pixel values outside epithlium to NaN


N = 20;
thresh = multithresh(IMgvec(isfinite(IMgvec)),N); % calculate N otsu threshROIds based on pixel values inside epithelium

disp('...')
disp('Optimist mask - Making a mask that defines what is definitely NOT a nucleus and should never be included.')

if manual_input_yes_no==1
    close all
    figure('units','normalized','outerposition',[0 0.2 1 0.6])
    for nn=1:N-10
        subplot(2,ceil((N-10)/2),nn)
        imshow(labeloverlay(IMg.*2,bwlabel(uint8(IMg>thresh(nn) & ROI==1),4),'Transparency',0.1));
        %title(['T(' num2str(nn)  ')=' num2str(round(thresh(nn),2))])
        title(['Image # ' num2str(nn) ]) 
        axis off
    end
    disp('Determine optimist threshold')
    prompt_opt = 'What thresholded image includes all nuclei and not more? Note: Nuclei should not have holes in them. (Enter image number from fig. here and press enter) : ';
    opt_index = input(prompt_opt);
elseif opt_index==0
    opt_index = 3;
end

disp('...')
disp('Pessimist mask, global - Making a mask that includes only the very brigtest pixels globally which should always be included. ')

if manual_input_yes_no==1
    close all
    % inspect how different thresholds perform and give manual input in command window
    figure('units','normalized','outerposition',[0 0.2 1 0.6])
    for nn=1:N-10
        subplot(2,ceil((N-10)/2),nn)
        imshow(labeloverlay(IMg.*2,bwlabel(uint8(IMg>thresh(nn) & ROI==1),4),'Transparency',0.1));
        %title(['T(' num2str(nn)  ')=' num2str(round(thresh(nn),2))])
        title(['Image # ' num2str(nn) ])
        axis off
    end
    disp('Determine pessimist threshold. ')
    prompt_pess = 'What thresholded image merges no nuclei? Very few merges are ok... (Enter image number from fig. here and press enter) : ';
    pess_index = input(prompt_pess);
elseif pess_index==0
    pess_index = 9;
end

optimist = IMg>thresh(opt_index); % this threshold should include everything
pessimist = IMg>thresh(pess_index); % this threshold should include only very brightest
pessimist = imfill(pessimist,'holes');

optimist(ROI==0)=0;
pessimist(ROI==0)=0;

regT = IMg-imgaussfilt(IMg,[15 15]);
regT = regT-min(regT(:));
regT = regT./(max(regT(:)));
regT(ROI==0)=NaN;

regTvec = regT(:);

N = 20;
thresh = multithresh(regTvec(isfinite(regTvec)),N); % calculate N otsu thresholds based on pixel values inside epithelium

disp('...')
disp('Regional threshold mask - Making a mask based on local thresholds decided by local ave. intensity')

if manual_input_yes_no==1
    close all
    % inspect how different thresholds perform and give manual input in command window
    figure('units','normalized','outerposition',[0 0.2 1 0.6])
    for nn=1:N-10
        subplot(2,ceil((N-10)/2),nn)
        imshow(labeloverlay(IMg.*2,bwlabel(uint8(imfill(regT>thresh(nn),'holes') & ROI==1 & optimist==1),4),'Transparency',0.1));
        %title(['T(' num2str(nn)  ')=' num2str(round(thresh(nn),2))])
        title(['Image # ' num2str(nn) ])
        
        axis off
    end
    disp('Fine tune regional thresholding.  ')
    prompt_regT = 'What thresholded image looks best?. I.e. separates most nuclei while losing the least? Smoother boundaries are better than seperation. (Enter image number from fig. here and press enter) : ';
    regT_index = input(prompt_regT);
elseif regT_index==0
    regT_index = 7;
end

regT = imfill(regT>thresh(regT_index),'holes');
regT(ROI==0)=0;
regT = bwareaopen(regT,approxPixelWidthOfCellNuclei,4);% remove small connected objects
regT(optimist==0)=0;% if its not in optimist mask remove

altPess = IMg-imgaussfilt(IMg,[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(25/50));
altPess = altPess-min(altPess(:));
altPess = altPess./(max(altPess(:)));
altPess(ROI==0)=NaN;

altPessvec = altPess(:);

N = 20;
thresh = multithresh(altPessvec(isfinite(altPessvec)),N); % calculate N otsu thresholds based on pixel values inside epithelium

disp('...')
disp('Pessimist mask, regional - Making a mask that includes only the very brigtest pixels locally, which should always be included. ')


if manual_input_yes_no==1
    % inspect how different thresholds perform and give manual input in command window
    close all
    figure('units','normalized','outerposition',[0 0.2 1 0.6])
    for nn=5:N-5-1
        subplot(2,ceil((N-10)/2),nn-5+1)
        imshow(labeloverlay(IMg.*2,bwlabel(uint8(imfill(altPess>thresh(nn),'holes') & ROI==1 & optimist==1),4),'Transparency',0.1));
        title(['Image # ' num2str(nn) ])        
        axis off
    end
    disp('Determine regional pessimist threshold: ')
    prompt_altPess ='What thresholded image merges no nuclei? Very few merges are ok... (Enter image number from fig. here and press enter) : ';
    altPess_index = input(prompt_altPess);
elseif altPess_index==0
    altPess_index = 10;
end

altPess = altPess>thresh(altPess_index);
altPess(ROI==0)=0;
altPess(optimist==0)=0;% if its not in optimist mask remove

pessimist = pessimist | altPess; % combine the two pessimists (then always include what the pessimist found)
pessimist = imfill(pessimist,'holes');

regT(pessimist==1)=1;% if pessimist found it then always include
regT = bwareaopen(regT,approxPixelWidthOfCellNuclei,4);% remove small connected objects
regT = imfill(regT,'holes');
regT(optimist==0)=0;

%figure
%imshow(labeloverlay(IMg,bwlabel(regT,4)) ); title(['regT: regional theshold, include nothing the optimist didnt find and include everything the pessimist did find'])

% Seperation of some merged nuclei via imopen (morphogical operation) and then smoothing of object boundaries to make later watershed perform better
% We use imopen more agressively on larger connected comp (larger object are more likely to be several nuclei in need of seperation)
disp('...')
disp('Doing imopen and smoothing on each connected component in the binary image regT seperately... this takes a while...')

L1 = bwlabel(regT);

TEMPim = uint8(zeros(size(L1)));

numCC = max(L1(:));

tic
for ii = 1:numCC% iterate through all connected components
    
    if mod(ii,100)==0
        disp([num2str(ii) ' out of ' num2str(numCC) ' objects imopened'])
    end
    
    y1o = approxPixelWidthOfCellNuclei*(1/50);% lower threshold for disk radius used for imopen on small objects
    y2o = approxPixelWidthOfCellNuclei*(7/50);% upper threshold for disk radius used for imopen on larger objects
    x1o = pi*(approxPixelWidthOfCellNuclei/4)^2; % everything below this object size thresh gets imopen with disk size y1o
    x2o = 2*x1o;% everything above this object size thresh gets imopen with disk size y2o
    ao = (y2o-y1o)/(x2o-x1o); % slop of liniar function deciding disk size for object with size between x1o and x2o
    bo = y1o - ao*x1o; % y axis intyersect of liniar function deciding disk size for object with size between x1o and x2o
    
    CurrentComIm = uint8(L1==ii); % image containing only current component/object
    currentx = sum(CurrentComIm(:)); % size of current component
    
    % find what disk size to use for imopen (floor(currenty))
    if currentx<=x1o
        currenty = x1o*ao + bo;
    elseif currentx>x1o && currentx<x2o
        currenty = currentx*ao + bo;
    elseif currentx>=x2o
        currenty = x2o*ao + bo;
    end
    
    % perform imopen with a strel that depends on the size of the connected
    % component in relation to the average pixel size of a nucleus (larger object are more likely to be several nuclei in need of seperation)
    TEMPim = TEMPim + imopen(CurrentComIm,strel('disk',floor(currenty)));% add new object one by one to empty image - whenever you get overlap between two you will get numbers greater than 1
    
end

TEMPim(TEMPim>1)=0;% remove regions that overlap after imopen (they have integers larger than 1)

% Smooth boundaries without merging connected comp

% some larger regions will have been split by imopen now so we relabel
% first
L1 = bwlabel(TEMPim,4);

TEMPim2 = uint8(zeros(size(L1)));% initialize empty temp image

numCC = max(L1(:));

for ii = 1:numCC% iterate through all objects
    
    if mod(ii,100)==0
        disp([num2str(ii) ' out of ' num2str(numCC) ' objects smoothed'])
    end
    CurrentComIm = uint8(L1==ii);% image containing only current component/object
    
    TEMPim2 = TEMPim2 + uint8(imgaussfilt(CurrentComIm,[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(2.5/50))>0.5);
    
end
disp('Done!')
toc

TEMPim2(TEMPim2>1)=0;% remove regions where smoothed regions overlap

regTcs = TEMPim2;% new name for imclosed smoothed regT!

% figure
% imshowpair(regTcs ,regT); title(['difference between regT and imclosed and smoothed regT'])
%
% figure
% imshow(labeloverlay(IM,bwlabel(regTcs,4)))

% Distance transform on regTcs and then watershed for seperation of nuclei that are still merged

disp('...')
disp('Now doing distance transform and watershed.')

inv_regTcs = -1.*double(regTcs)+1;

D = bwdist(inv_regTcs);% distance to nucleus edge

Dg = imgaussfilt(D,[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(3/50));% gauss filter dist image to reduce the number of 'false' regional max

diskRadius1 = round(approxPixelWidthOfCellNuclei.*(1/50));

rm = imregionalmax(Dg);% locate regional max pixels in smoothed distance transfrom image
rm = imdilate(rm,strel('disk',diskRadius1));% dilate regional max a bit
diskRadius2 = round(approxPixelWidthOfCellNuclei.*(3/50));
rm = imclose(rm,strel('disk',diskRadius2));% imclose rm - merges regional max which are very close together (reducing over segmentation)
rm = imdilate(rm,strel('disk',diskRadius1));% imclose rm - merges regional max which are very close together (reducing over segmentation)
rm(regTcs==0)=0; % remove pixels if they where dilated outside regTcs

w = watershed(imimposemin(-1.*D, rm | inv_regTcs));% impose min on -1.*D at rm and inverse regTcs and do watershed on -1.*D

% process watershed result
FC = w>0; % isolate watershed lines
FC(regTcs==0)=0;

range = [2.*approxPixelWidthOfCellNuclei 2*pi*(approxPixelWidthOfCellNuclei/2)^2];% exclude object smaller than 2.*approxPixelWidthOfCellNuclei and larger than 2 times ave nucleus area
FC = bwareafilt(FC,range,4);
%FC = bwareaopen(FC,2.*approxPixelWidthOfCellNuclei,4);

BW = logical(FC);

% Comment out if you dont want to see final result plotted!
figure
imshow(labeloverlay((IM./max(IM(:))),bwlabel(BW,4))); title('Final nuclei segmentation. If result is oversegmented then try again with a slighty larger approx. pixel diameter of nucleus. Else set manual-input-yes-no==1 and answer prompts in command line')

end

