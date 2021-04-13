close all
clear all

path = pwd;  % pwd returns the path to the current folder.
subfolderLSM = '/raw_images_lsm_files';
subfolderMAT = '/outlines_mat_files';

list_lsm = dir([path filesep subfolderLSM filesep '*.lsm']);
list_mat = dir([path filesep subfolderMAT filesep '*OL.mat']);

ff = 7;%% Run with manual input for just one image file in folder (raw_images_lsm_files) with 8 test images

lsm_image_string   = [path filesep subfolderLSM filesep list_lsm(ff).name];
mat_outLine_string = [path filesep subfolderMAT filesep list_mat(ff).name];

% ------------------------------------------------------------------------    
% Load image file and handdrawn outline    

load(mat_outLine_string);% load hand drawn epithelial outline saved in .mat file

t = Tiff(lsm_image_string); % create tiff object
I = read(t); % read in image file (Tiff() reader works on .lsm zeiss files)
close(t);% close tiff object

C1 = double(I(:,:,1));% make double
C2 = double(I(:,:,2));
C3 = double(I(:,:,3));
C4 = double(I(:,:,4));

% Make normalized images for viewing  (normalize with respect to max
% value)
C1norm = C1./max(C1(:));
C2norm = C2./max(C2(:));
C3norm = C3./max(C3(:));
C4norm = C4./max(C4(:));

% make normalized rgb color image for quick viewing
COL = uint8(zeros(size(I,1),size(I,2),3));
COL(:,:,1) = uint8(C3norm.*255) + uint8(C1norm.*255).*5;% red Sox9 (C3)      , white CPA (C1)
COL(:,:,2) = uint8(C2norm.*255) + uint8(C1norm.*255).*5;% green p120ctn (C2) , white CPA (C1)
COL(:,:,3) = uint8(C4norm.*255) + uint8(C1norm.*255).*5;% blue DAPI (C4)     , white CPA (C1)

OL = imfill(OL,'holes');

% make three layer segmentation edge image for inclusion in color image (COL_show)
% when shown
OL3edge = uint8(zeros(size(OL,1),size(OL,2),3)); OL3edge(:,:,1) = bwperim(OL); OL3edge(:,:,2) = bwperim(OL); OL3edge(:,:,3) = bwperim(OL);

COL_show = COL; COL_show(imdilate(OL3edge,strel('disk',3))==1)=2^8-1;


% ------------------------------------------------------------------------    
% Tip score calculation from outline OL    

figure('units','normalized','outerposition',[0 0.2 1 0.6]);
imshow(COL_show,[]); title('Red: Ecad. Blue: DAPI. Green: p120ctn. White: CPA  - Measure the average length of tip structures in pixels (needed as manual input)')
h = imdistline;
fcn = makeConstrainToRectFcn('imline', get(gca,'XLim'),get(gca,'YLim'));
setDragConstraintFcn(h,fcn);
disp('')
prompt_tip = 'What is the average diameter/length of tip structures in pixels? (~300 for test images) Enter integer number here and press enter) : ';
approxPixelWidthOfTipStructure = round(input(prompt_tip));

plot_yes_no = 1; % if set to 1 you will see a plot of the tip score and its components (seeing this will help decide if the chosen approxPixelWidthOfTipStructure is too large/small)

% ##################################################################
TSC = tipScoreIm(OL,approxPixelWidthOfTipStructure,plot_yes_no); % tipscore function is calles and returns an image with tipscores
% ##################################################################


% ------------------------------------------------------------------------    
% DAPI image segmentation

% find ave length of nucleus - this is mandatory input for the nucleus segmentation function segmentDAPIimage()

figure('units','normalized','outerposition',[0 0.2 1 0.6]);
imshow(COL_show,[]); title('Red: Ecad. Blue: DAPI. Green: p120ctn. White: CPA  - Measure the average diameters of nuclei in pixels (needed as manual input)')
h = imdistline;
fcn = makeConstrainToRectFcn('imline', get(gca,'XLim'),get(gca,'YLim'));
setDragConstraintFcn(h,fcn);
disp('')
prompt_nuc = 'What is the average diameter of nuclei in pixels? (15-20 for test images) Enter integer number here and press enter) : ';
approxPixelWidthOfCellNuclei = round(2.*input(prompt_nuc));
manual_input_yes_no = 1; 
IM = C4;% DAPI channel!
ROI = OL; % mask of eg. epithelial outline. If you want to segment nuclei everywhere then input image with all zeros.
answers = [NaN NaN NaN NaN];

% ##################################################################
NUCLEI_BW = segmentDAPIimage(IM,ROI,approxPixelWidthOfCellNuclei,manual_input_yes_no,answers);% Segment nuclei in DAPI image, returns bimary image
% ##################################################################

figure
imshow(labeloverlay((IM./max(IM(:))),bwlabel(NUCLEI_BW,4))); title('Final nuclei segmentation. If result is oversegmented then try again with a slighty larger approx. pixel diameter of nucleus. Else set manual-input-yes-no==1 and answer prompts in command line')


% ---------------------------------------------------------------------------------------------------------------------
% Make Tcells table with image intensities and tip scores per cell for image using the returnTableWithCellInt() function

approxPixelWidthOfCellNuclei = 20; % this desides the level of gauss smoothing of the image

TcellsC1 = returnTableWithCellInt(NUCLEI_BW,C1,approxPixelWidthOfCellNuclei);
TcellsC1.Properties.VariableNames{'IMintensity'} = 'CPAintensity';

TcellsC3 = returnTableWithCellInt(NUCLEI_BW,C3,approxPixelWidthOfCellNuclei);
TcellsC3.Properties.VariableNames{'IMintensity'} = 'Sox9intensity';

approxPixelWidthOfCellNuclei = 1;% we dont wnt to gauss smooth tipscore image so set to 1
TcellsTSC = returnTableWithCellInt(NUCLEI_BW,TSC,approxPixelWidthOfCellNuclei);

TcellsTSC.Properties.VariableNames{'IMintensity'} = 'tipScore';

Tcells = [TcellsC1,TcellsC3(:,4),TcellsTSC(:,4)];% merge tables

% ---------------------------------------------------------------------------------------------------------------------
% normalize intensities with respect to median intensity in a group of
% cells with similar tip scores (for comparison with intensities in
% other images).
Tcells.CPAintensityNorm  = Tcells.CPAintensity./nanmedian(Tcells.CPAintensity(Tcells.tipScore>-0.75 & Tcells.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
Tcells.Sox9intensityNorm = Tcells.Sox9intensity./nanmedian(Tcells.Sox9intensity(Tcells.tipScore>-0.75 & Tcells.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)

Tcells(1:3,:)% show first three rows of final table in command line

% ---------------------------------------------------------------------------------------------------------------------
% Make TcellPairs table with image intensities and tip scores for image using the returnTableWithCellPairInt() function

approxPixelWidthOfCellNuclei = 35; % sets threshold for which cells are considered direct neighbors (larger number means more neiborg pairs and increased chance of line scab crossing two cell membrane)
plot_yes_no = 1; % WARNING!! - drawing a lot of lines on an image takes time so set this to 0 when you dont need to inspect where the line scans are done!!!

TcellPairsp120ctn  = returnTableWithCellPairInt(NUCLEI_BW,C2,approxPixelWidthOfCellNuclei,plot_yes_no);
TcellPairsp120ctn.Properties.VariableNames{'lineScanMax'} = 'p120ctnMaxInt';% change genric name to specific marker name
TcellPairsp120ctn.Properties.VariableNames{'lineScanMed'} = 'p120ctnMedInt';% change genric name to specific marker name
TcellPairsp120ctn.Properties.VariableNames{'lineScan'}    = 'p120ctnLineScan';

TcellPairsTSC = returnTableWithCellPairInt(NUCLEI_BW,TSC,approxPixelWidthOfCellNuclei,0);
TcellPairsTSC.Properties.VariableNames{'lineScanMed'} = 'tipScore';% we only used the median for the tip score of the cell pairs

TcellPairs = [TcellPairsp120ctn,TcellPairsTSC(:,9)];% merge tables

TcellPairs.p120ctnMaxIntNorm = TcellPairs.p120ctnMaxInt./nanmedian(TcellPairs.p120ctnMaxInt(TcellPairs.tipScore>-0.75 & TcellPairs.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
TcellPairs.p120ctnMedIntNorm = TcellPairs.p120ctnMedInt./nanmedian(TcellPairs.p120ctnMaxInt(TcellPairs.tipScore>-0.75 & TcellPairs.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
  
% ---------------------------------------------------------------------------------------------------------------------
% Make coarse grained tip score for grouping (the lower factor is - the more groups)
factor = 3;
TcellPairs.tipScoreGroup =  round(TcellPairs.tipScore./factor,1)*factor;
Tcells.tipScoreGroup     =  round(Tcells.tipScore    ./factor,1)*factor;

tipscoregroups = unique(TcellPairs.tipScoreGroup);
tipscoregroups = tipscoregroups(isfinite(tipscoregroups));

% ---------------------------------------------------------------------------------------------------------------------
% the dist. tails that have very low or very high tip scores contain few
% meassurements - pool them
index1 = 6;
TcellPairs.tipScoreGroup(TcellPairs.tipScoreGroup<tipscoregroups(index1)) = tipscoregroups(index1);
Tcells.tipScoreGroup(Tcells.tipScoreGroup<tipscoregroups(index1)) = tipscoregroups(index1);
index2 = length(tipscoregroups)-1;
TcellPairs.tipScoreGroup(TcellPairs.tipScoreGroup>tipscoregroups(index2)) = tipscoregroups(index2);
Tcells.tipScoreGroup(Tcells.tipScoreGroup>tipscoregroups(index2)) = tipscoregroups(index2);

TcellPairs(1:3,:)% write out first 3 rows of table to command line


% ---------------------------------------------------------------------------------------------------------------------
% Plot Sox9 intensity per cell as a function of tipscore. And compare groups with Kruskal Wallis
figure
boxplot(Tcells.Sox9intensityNorm,Tcells.tipScoreGroup,'Notch','on','jitter',.1,'symbol','.'); hold on
xlabel('Tip score'); set(gca,'xticklabel',{'TS below -0.6','TS between -0.6 and -0.2','TS between -0.2 and 0.2','TS above 0.2' })

ylabel('Normalized  Sox9 intensity per cell')

[pSox9,tblSox9,statsSox9] = kruskalwallis(Tcells.Sox9intensityNorm,Tcells.tipScoreGroup);
[cSox9,mSox9,hSox9,gnamesSox9] = multcompare(statsSox9);

% ---------------------------------------------------------------------------------------------------------------------
% Plot p120ctn intensity as a function of (coarse grained) tip score. And compare groups with Kruskal Wallis
figure 
boxplot(TcellPairs.p120ctnMaxIntNorm,TcellPairs.tipScoreGroup,'Notch','on','jitter',.1,'symbol','.'); hold on
xlabel('Tip score');  set(gca,'xticklabel',{'TS below -0.6','TS between -0.6 and -0.2','TS between -0.2 and 0.2','TS above 0.2' })

ylabel('Normalized max p120ctn int. per cell-pair linescan')

[pP120,tblP120,statsP120] = kruskalwallis(TcellPairs.p120ctnMaxIntNorm,TcellPairs.tipScoreGroup)
[cP120,mP120,hP120,gnamesP120] = multcompare(statsP120);

