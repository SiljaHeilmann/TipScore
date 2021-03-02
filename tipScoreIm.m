function [tipScoreIm] = tipScoreIm(BranchedMask,approxPixelWidthOfTipStructure,plot_yes_no)
% tipScoreIm takes a BW mask of a branched structure approx width/length of tip structures
% in the tree and returns an image (of type double) where each pixel
% contains a value that quantifies how 'tip like' its
% location is, (values >0 will usually be in quite tip-like areas). (Tip score 
% is usually in the range [-2,1], but number depend on the input value 'approxPixelWidthOfTipStructure' 
% and the specific structure given by 'BranchedMask'). Pixel values outside the mask is set to NaN. 
% The tip score function was developed to quantify/correlate well with
% what biologist judge by eye as 'tip' (and anti tip = trunk) in the
% branched structure of a pancreas (but could potentially be used to
% quantify location in other branched structures too).
%
% INPUTS:
% BranchedMask: a binary mask (logical, double, uint8, or uint16) defining the branched
% structure. 1's everywhere inside structure, 0's oustside. THis mask could be
% based on a handdrawn outline or based on a segmentation.
% approxPixelWidthOfTipStructure: the lengthscale/width/length measured in pixels of what you
% judge to be a 'tip structure' in your image.
% plot_yes_no: if set to 1 you will see a plot of the final tipscore image
% together with a plot of the components used for calculate it. (Ei.
% sidebranch, midbranch and convex hull images). Seeing this can help you
% fine tune the 'approxPixelWidthOfTipStructure' to an appropriate value.
%
% OUTPUT:
% tipScoreIm: image (of type double) where each pixel
% contains a value that quantifies how 'tip like' its
% location is, (values >0 will be in quite tip-like areas). (Tip score 
% is usually in the range [-2,1], but number depend on the input value 'approxPixelWidthOfTipStructure' 
% and the specific structure given by 'BranchedMask'). Pixel values outside the mask is set to NaN.


disp('.')
disp('Starting...')

plot_yes = plot_yes_no;

disp('...')
disp('Skeletonize branching structure.')

BranchedMask = logical(BranchedMask);

OLf = medfilt2(imdilate(uint8(imfill(BranchedMask,'holes')),strel('disk',1)),[3 3]);% dilate one pixel and median filter to reduce spurs in skeleton, fill small holes (connected structure needs to be hole free!)

SKspur =  logical(bwmorph(logical(OLf),'thin',Inf));% Thin structure to single pixel width skeleton
MustIncludePoints = bwmorph(SKspur,'shrink',Inf);% Thin structure to single points (for safe keeping to not loose small parts when skeleton is reduced later)
SK = bwskel(SKspur,'MinBranchLength',round(approxPixelWidthOfTipStructure.*(1/4))); % remove very small spurs/sidebranches
SK = bwskel(SK);% skeletonize again to avoid that some inner points get labelled as branch points wrongly

SK(MustIncludePoints==1)=1;% include these points if they where lost when removing small sidebranches (to avoid loosing an entire connected comp that was small...)

BP = bwmorph(SK,'branchpoint');% find branchpoints

SKred = SK; % define new skeleton image for step wise reduction

for ii = 1:approxPixelWidthOfTipStructure % we define length of side branch as the length of an ave. tip structure
    
    EP_temp = bwmorph(SKred,'endpoint');% find current end points in reduced skeleton SKred
    EP_temp(BP==1)=0; % dont an endpoint if it is branchpoint!
    SKred = SKred-EP_temp;% remove endpoints from skeleton (thereby shortening branches)
    
end

disp('...')
disp('Finding side branches and midbranch.')

MB = SKred;% Reduced skeleton is midbranch
SB = SK-SKred; % The part that was removed from skeleton is sidebranch

LOLf = bwlabel(OLf);% label individual connected components in OLf such that we may find the convex hull of each individual part
CH = zeros(size(OLf)); % initialize image for convex hull
CHdist = zeros(size(OLf)); % initialize image for distance to convex hull
MBdist = zeros(size(OLf)); % initialize image for distance to midbranch
SBdist = zeros(size(OLf)); % initialize image for distance to sidebranch

disp('...')
disp('Determine convex hull of each connected component.')


for ii=1:max(LOLf(:))% iterate through the connected components
    TEMP1 = zeros(size(OLf));% initialize temporary image
    TEMP2 = zeros(size(OLf));% initialize temporary image
    
    r = regionprops(LOLf,'Image','ConvexImage','BoundingBox');% get Image, ConvexImage and BoundingBox of each connected component in OLf
    BB = ceil(r(ii).BoundingBox);% round of to integer
    
    TEMP2(BB(2):BB(2)+BB(4)-1,BB(1):BB(1)+BB(3)-1)=bwperim(r(ii).ConvexImage); % insert edge of convex hull in larger empty temp image
    
    CH = CH+TEMP2;% add current convex hull edge to CH
    
    IM = bwdist(-1.*r(ii).ConvexImage+1);% find distance to convex hull edge
    IM(r(ii).Image==0) = 0;% set values outside epithelium to 0
    TEMP1(BB(2):BB(2)+BB(4)-1,BB(1):BB(1)+BB(3)-1)=IM; % insert dist to convex hull in larger empty temp image
    CHdist = CHdist+TEMP1; % add current dist to convex hull to CHdist
    
    currentMB = MB;
    currentMB(LOLf~=ii)=0;% make image with only the mid branches from the current connected comp
    currentSB = SB;
    currentSB(LOLf~=ii)=0;% make image with only the side branches from the current connected comp
    
    currentMBdist = bwdist(currentMB);% make image with only the distance to the mid branches from the current connected comp
    currentMBdist(~isfinite(currentMBdist))=0;% do not include Inf values (if current comp doesnt have a midbranch we get Inf)
    currentMBdist(LOLf~=ii)=0;% set values outside current connected comp to 0
    
    currentSBdist = bwdist(currentSB);% make image with only the distance to the side branches from the current connected comp
    currentSBdist(~isfinite(currentSBdist))=0;% do not include Inf values (if current comp doesnt have a side branch we get Inf)
    currentSBdist(LOLf~=ii)=0;% set values outside current connected comp to 0
    
    MBdist = MBdist + currentMBdist;% add current dist to mid branch to MBdist
    SBdist = SBdist + currentSBdist;% add current dist to side branch to SBdist
    
end

disp('...')
disp('Calculating final tip score.')


% distance values are normalised with length scale of tip structure
CHsq_norm = sqrt(CHdist./approxPixelWidthOfTipStructure);% we take the sqrt of the distance to convex hull
MB_norm = MBdist./approxPixelWidthOfTipStructure;
SB_norm = SBdist./approxPixelWidthOfTipStructure;

SB_norm(OLf==0) = NaN;
MB_norm(OLf==0) = NaN;
CHsq_norm(OLf==0) = NaN;

% calculate tic score image
TSC = MB_norm - SB_norm - CHsq_norm;

MAX = 2;% puts an upper and lower limit on tip score values (comment out if needed)
TSC(TSC>MAX)=MAX;
TSC(TSC<-MAX)=-MAX;

TSC(OLf==0)=NaN;% set values outside epithelial mask to NaN

if plot_yes==1
    figure('units','normalized','outerposition',[0 0.2 1 0.6])%[0 0 1 1]
    subplot(241)
    imshowpair(imdilate(bwperim(CH>0),strel('disk',5)),imdilate(bwperim(OLf>0),strel('disk',5))); axis off; title('CH')
    subplot(242)
    imshowpair(imdilate(MB,strel('disk',5)),imdilate(bwperim(OLf>0),strel('disk',5))); title('MB'); axis off;
    subplot(243)
    imshowpair(imdilate(SB,strel('disk',5)),imdilate(bwperim(OLf>0),strel('disk',5))); title('SB'); axis off;
%    subplot(244)
%    imshow(COL_show,[]);title('Image. Red: Ecad. Blue: DAPI. Green: p120ctn. White: CPA'); axis off;
    
    subplot(245)
    imshow(CHsq_norm,[]); title('sqrt(dist to CH)'); axis off;
    subplot(246)
    imshow(MB_norm,[]); title('dist to MB'); axis off;
    subplot(247)
    imshow(SB_norm,[]); title('dist to SB'); axis off;
    subplot(248)
    imshow(TSC,[]); title('Tipscore = dist(MB) - dist(SB) - sqrt[dist(CH)]'); axis off;
    colormap jet
    
end

tipScoreIm = TSC;



end

