close all

path = pwd;  % pwd returns the path to the current folder.
subfolderTIF = 'raw_images_tif_files';
subfolderMAT  = 'outlines_mat_files';

list_tif = dir([path filesep subfolderTIF filesep '*.tif']);
list_mat  = dir([path filesep subfolderMAT filesep '*OLnew.mat']);

clear DS
DS(1).ImageFileName = "";% initialize data structure for storing images, masks and extracted data of several images

% input for tipScoreIm() function
plot_yes_no = 0;
approxPixelWidthOfTipStructure = 300;

% input for segmentDAPIimage() function
approxPixelWidthOfCellNuclei = 20;% resolution of all images is the same so we can use th esame length for all images in the folder
manual_input_yes_no = 0;
 
   answersMatrix = ... % matrix containing appropriate answers for the four promts given when the segmentDAPIimage() function is run with 'manual_input_yes_no=1'
    [4 7 6 9
     4 7 6 9
     4 7 6 9
     4 7 6 9
     4 7 6 9
     4 7 6 9
     4 7 6 9];
 
 
for ff = 1:length(list_tif)
    
    tif_image_string   = [path filesep subfolderTIF filesep list_tif(ff).name];
    mat_outLine_string  = [path filesep subfolderMAT  filesep list_mat(ff).name];
    
    load(mat_outLine_string);% load hand drawn epithelial outline saved in .mat file
    
    OL = OLnew;% OL is called OLnew in .mat file
    clear OLnew

    C1 = double(imread(tif_image_string,1));% make double
    C2 = double(imread(tif_image_string,2));
    C3 = double(imread(tif_image_string,3));
    C4 = double(imread(tif_image_string,4));
    
    % Make normalized images for viewing  (normalize with respect to 99
    % percentile value)
    C1norm = C1./prctile(C1(:),99.0);
    C2norm = C2./prctile(C2(:),99.0);
    C3norm = C3./prctile(C3(:),99.0);
    C4norm = C4./prctile(C4(:),99.0);

    % make normalized rgb color image for quick viewing
    COL = uint8(zeros(size(C1,1),size(C1,2),3));

    COL(:,:,1) = uint8(C3norm.*255) + uint8(C1norm.*255);% red Sox9 (C3)      , white CPA (C1)
    COL(:,:,2) = uint8(C2norm.*255) + uint8(C1norm.*255);% green p120ctn (C2) , white CPA (C1)
    COL(:,:,3) = uint8(C4norm.*255) + uint8(C1norm.*255);% blue DAPI (C4)     , white CPA (C1)

    OL = imfill(OL,'holes');
    
    % make three layer segmentation edge image for inclusion in color image (COL_show)
    % when shown
    OL3edge = uint8(zeros(size(OL,1),size(OL,2),3)); OL3edge(:,:,1) = bwperim(OL); OL3edge(:,:,2) = bwperim(OL); OL3edge(:,:,3) = bwperim(OL);
    
    COL_show = COL;
    COL_show(imdilate(OL3edge,strel('disk',3))==1)=2^8-1;
    
    imageName = list_tif(ff).name;
    imageFolderForOutput = [imageName(1:end-4) '_Image#' num2str(ff) ];
    
    % make new folder for output
    command = ['mkdir  -p ' '"' pwd filesep imageFolderForOutput '"']
    
    status = system(command,'-echo'); % echo command line in matlab command window
    
    % save outline OL and each channel (C1,C2,C3 and C4 plus COL_show) as mat files in new folder for quicker access later
    save([imageFolderForOutput filesep 'OL.mat'],'OL')
    save([imageFolderForOutput filesep 'C1.mat'],'C1')
    save([imageFolderForOutput filesep 'C2.mat'],'C2')
    save([imageFolderForOutput filesep 'C3.mat'],'C3')
    save([imageFolderForOutput filesep 'C4.mat'],'C4')
    save([imageFolderForOutput filesep 'COL_show.mat'],'COL_show')
    
    disp(['Now processing image number ' num2str(ff) ' out of ' num2str(length(list_tif))])

% ------------------------------------------------------------------------    
% Tip score calculation from outline OL    
    
    % ##################################################################
    TSC = tipScoreIm(OL,approxPixelWidthOfTipStructure,plot_yes_no);
    % ##################################################################
    
    save([imageFolderForOutput filesep 'TSC.mat'],'TSC')
    
% ------------------------------------------------------------------------    
% DAPI image segmentation

    IM = C4;% Here channal  4 is DAPI image
    
    ROI = OL; % mask of eg. epithelial outline. If you want to segment nuclei everywhere then input ROI = image with all zeros.
    
    answers = answersMatrix(ff,:);% can be used to give unique answer per image
    
    % ##################################################################
    NUCLEI_BW = segmentDAPIimage(IM,ROI,approxPixelWidthOfCellNuclei,manual_input_yes_no,answers);
    % ##################################################################
    
    save([imageFolderForOutput filesep 'NUCLEI_BW.mat'],'NUCLEI_BW')
    
% ------------------------------------------------------------------------    
% Store everything in data structure DS
        
    DS(ff).OrgPath       = string(path);
    DS(ff).OrgSubfolder  = string(subfolderTIF);
    DS(ff).ImageFileName = string(imageName);
    DS(ff).FolderName    = string(imageFolderForOutput);
    DS(ff).COL_show      = COL_show;
    DS(ff).TSC           = TSC;
    DS(ff).NUCLEI_BW     = NUCLEI_BW;
    DS(ff).OL            = OL;
    DS(ff).C1            = C1;
    DS(ff).C2            = C2;
    DS(ff).C3            = C3;
    DS(ff).C4            = C4;
    
end

save([pwd filesep 'DS.mat'],'DS','-v7.3')% save data structure in current dir

%load([pwd filesep 'DS.mat'])

% ------------------------------------------------------------------------    
% Plot/show tip score images for each image in the test set in one figure
% close all
figure
for ff=1:length(DS)
    subplot(ceil(sqrt(length(DS))),round(sqrt(length(DS))),ff)
    %imshowpair(DS(ff).C4./max(DS(ff).C4(:)),DS(ff).NUCLEI_BW);
    imshow(DS(ff).TSC,[-2 1]);
    colorbar
    colormap jet
    axis off
    title(string([num2str(ff) 'Filename: '  DS(ff).ImageFileName]));
end


% ------------------------------------------------------------------------    
% Make one Tcells table with image intensities and tip scores for each image using the returnTableWithCellInt() function
% adds a column with normalized intensities. Normalization is done with
% respect to intensities of a group of cells that have similar tipscores

for ff=1:length(DS)

    disp([ 'Image # ' num2str(ff) ' out of ' num2str(length(DS)) ]);
    
    approxPixelWidthOfCellNuclei = 27;
    
    TcellsC1 = returnTableWithCellInt(DS(ff).NUCLEI_BW,DS(ff).C1,approxPixelWidthOfCellNuclei);
    TcellsC1.Properties.VariableNames{'IMintensity'} = 'CPAintensity';
    
    TcellsC3 = returnTableWithCellInt(DS(ff).NUCLEI_BW,DS(ff).C3,approxPixelWidthOfCellNuclei);
    TcellsC3.Properties.VariableNames{'IMintensity'} = 'Sox9intensity';
    
    approxPixelWidthOfCellNuclei = 1;% we dont wnt to gauss smooth tipscore image so set to 1
    TcellsTSC = returnTableWithCellInt(DS(ff).NUCLEI_BW,DS(ff).TSC,approxPixelWidthOfCellNuclei);

    TcellsTSC.Properties.VariableNames{'IMintensity'} = 'tipScore';

    Tcells = [TcellsC1,TcellsC3(:,4),TcellsTSC(:,4)];% merge tables
    
    % normalize intensities with respect to median intensity in a group of
    % cells with similar tip scores (for comparison with intensities in
    % other images).
    Tcells.CPAintensityNorm  = Tcells.CPAintensity./nanmedian(Tcells.CPAintensity(Tcells.tipScore>-0.75 & Tcells.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
    Tcells.Sox9intensityNorm = Tcells.Sox9intensity./nanmedian(Tcells.Sox9intensity(Tcells.tipScore>-0.75 & Tcells.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)

    Tcells(1:3,:)
    
    DS(ff).Tcells = Tcells;% put table in data structure
    
end
disp('Done! Saving...')

save([pwd filesep 'DS.mat'],'DS','-v7.3')


% ------------------------------------------------------------------------    
% Make one TcellPairs table with image intensities and tip scores for each image
% adds a column with normalized intensities. Normalization is done with
% respect to intensities of a group of cell-cell pairs that have similar tipscores

plot_yes_no = 0;

for ff=1:length(DS)

    disp([ 'Image # ' num2str(ff) ' out of ' num2str(length(DS)) ]);
 
    approxPixelWidthOfCellNuclei = 35; % sets threshold for which cells are considered direct neighbors (larger number means more neiborg pairs and increased chance of line scab crossing two cell membrane)
 
    TcellPairsp120ctn  = returnTableWithCellPairInt(DS(ff).NUCLEI_BW,DS(ff).C2,approxPixelWidthOfCellNuclei,plot_yes_no);
    TcellPairsp120ctn.Properties.VariableNames{'lineScanMax'} = 'p120ctnMaxInt';% change genric name to specific marker name
    TcellPairsp120ctn.Properties.VariableNames{'lineScanMed'} = 'p120ctnMedInt';% change genric name to specific marker name
    TcellPairsp120ctn.Properties.VariableNames{'lineScan'}    = 'p120ctnLineScan';

    TcellPairsTSC = returnTableWithCellPairInt(DS(ff).NUCLEI_BW,DS(ff).TSC,approxPixelWidthOfCellNuclei,plot_yes_no);
    TcellPairsTSC.Properties.VariableNames{'lineScanMed'} = 'tipScore';% we only used the median for the tip score of the cell pairs
    
    TcellPairs = [TcellPairsp120ctn,TcellPairsTSC(:,9)];% merge tables

    TcellPairs.p120ctnMaxIntNorm = TcellPairs.p120ctnMaxInt./nanmedian(TcellPairs.p120ctnMaxInt(TcellPairs.tipScore>-0.75 & TcellPairs.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
    TcellPairs.p120ctnMedIntNorm = TcellPairs.p120ctnMedInt./nanmedian(TcellPairs.p120ctnMaxInt(TcellPairs.tipScore>-0.75 & TcellPairs.tipScore<-0.25));% Normalize intensity with respect median signal for group of cells with similar tipscore (tipScore=[-0.75;-0.25]is group that usually have a lot of cells...)
       
    TcellPairs(1:3,:)
    
    DS(ff).TcellPairs = TcellPairs;% put table in data structure

end

disp('Done! Saving...')

save([pwd filesep 'DS.mat'],'DS','-v7.3')

% ------------------------------------------------------------------------    
% Merge Tcells tables and TcellPairs tables for each image in two large tables: AllCells and AllCellPairs

AllCells = table;
AllCellPairs = table;

numCells = 0;
numCellPairs = 0;

ffcells = [];
ffcellPairs = [];


for ff = 1:length(DS)
    ff
    numCells     =  size(DS(ff).Tcells,1);
    numCellPairs =  size(DS(ff).TcellPairs,1);
    
    ffcells     = [ffcells,ones(1,numCells).*ff];
    ffcellPairs = [ffcellPairs,ones(1,numCellPairs).*ff];
    
    AllCells = [AllCells;DS(ff).Tcells];
    AllCellPairs = [AllCellPairs;DS(ff).TcellPairs];
    
end

% ------------------------------------------------------------------------    
% Make coarse grained tip score for grouping (the lower factor - the more groups)
factor = 3;
AllCellPairs.tipScoreGroup =  round(AllCellPairs.tipScore./factor,1)*factor;
AllCells.tipScoreGroup     =  round(AllCells.tipScore    ./factor,1)*factor;

tipscoregroups = unique(AllCellPairs.tipScoreGroup);
tipscoregroups = tipscoregroups(isfinite(tipscoregroups));

index1 = 5; 
tipscoregroups(index1)
AllCellPairs.tipScoreGroup(AllCellPairs.tipScoreGroup<tipscoregroups(index1)) = tipscoregroups(index1);% -1.5;
AllCells.tipScoreGroup(AllCells.tipScoreGroup<tipscoregroups(index1)) = tipscoregroups(index1);
 
index2 = length(tipscoregroups)-1;

AllCellPairs.tipScoreGroup(AllCellPairs.tipScoreGroup>tipscoregroups(index2)) = tipscoregroups(index2);
AllCells.tipScoreGroup(AllCells.tipScoreGroup>tipscoregroups(index2)) = tipscoregroups(index2);

AllCells.ff = ffcells';
AllCellPairs.ff = ffcellPairs';

save([pwd filesep 'AllCells.mat'],'AllCells')
save([pwd filesep 'AllCellPairs.mat'],'AllCellPairs')


% ------------------------------------------------------------------------    
%load([pwd filesep 'AllCells.mat'])
%load([pwd filesep 'AllCellPairs.mat'])
%load([pwd filesep 'DS.mat'])

% ------------------------------------------------------------------------    
% Plot Sox9 intensity as a function of (coarse grained) tip score. And compare groups with Kruskal Wallis

% figure
% for ff=1:length(DS)
%     ff
%     subplot(2,4,ff)
%     includeff = AllCells.ff==ff;
%     boxplot(AllCells.Sox9intensityNorm(includeff),AllCells.tipScoreGroup(includeff),'Notch','on','jitter',.1,'symbol','.'); hold on
%     axis([0 7 0 max(AllCells.Sox9intensityNorm)])
% %    xlabel('Tip score');  set(gca,'xticklabel',{'-0.9: (TS below -0.75)','-0.6: (TS between -0.75 and -0.45)','-0.3: (TS between -0.45) and -0.15','0: (TS between -0.15 and 0.15)','0.3: (TS between 0.15 and 0.45)','0.6: (TS above 0.45)' })
%     title(['Image #' num2str(ff)])
% end

figure
boxplot(AllCells.Sox9intensityNorm,AllCells.tipScoreGroup,'Notch','on','jitter',.1,'symbol','.'); hold on
xlabel('Tip score','Fontsize',15);   set(gca,'xticklabel',{'-0.9: (TS below -0.75)','-0.6: (TS between -0.75 and -0.45)','-0.3: (TS between -0.45) and -0.15','0: (TS between -0.15 and 0.15)','0.3: (TS between 0.15 and 0.45)','0.6: (TS above 0.45)' })


set(findobj(gca,'type','line'),'linew',2)

ylabel('Normalized  Sox9 intensity per cell','Fontsize',15)

[pSox9,tblSox9,statsSox9] = kruskalwallis(AllCells.Sox9intensityNorm,AllCells.tipScoreGroup);
[cSox9,mSox9,hSox9,gnamesSox9] = multcompare(statsSox9,'Alpha',0.05);


% ------------------------------------------------------------------------    
% Plot p120ctn intensity as a function of (coarse grained) tip score. And compare groups with Kruskal Wallis

% figure
% for ff=1:length(DS)
%     ff
%     subplot(2,4,ff)
%     includeff = AllCellPairs.ff==ff;
%     boxplot(AllCellPairs.p120ctnMaxIntNorm(includeff),AllCellPairs.tipScoreGroup(includeff),'Notch','on','jitter',.1,'symbol','.'); hold on
%     axis([0 7 0 max(AllCellPairs.p120ctnMaxIntNorm)])
%     xlabel('Tip score');  % set(gca,'xticklabel',{'-0.9: (TS below -0.75)','-0.6: (TS between -0.75 and -0.45)','-0.3: (TS between -0.45) and -0.15','0: (TS between -0.15 and 0.15)','0.3: (TS between 0.15 and 0.45)','0.6: (TS above 0.45)' })
%     title(['Image #' num2str(ff)])
% end

figure
boxplot(AllCellPairs.p120ctnMaxIntNorm,AllCellPairs.tipScoreGroup,'Notch','on','jitter',.1,'symbol','.'); hold on
xlabel('Tip score','Fontsize',15);  set(gca,'xticklabel',{'-0.9: (TS below -0.75)','-0.6: (TS between -0.75 and -0.45)','-0.3: (TS between -0.45) and -0.15','0: (TS between -0.15 and 0.15)','0.3: (TS between 0.15 and 0.45)','0.6: (TS above 0.45)' })
set(findobj(gca,'type','line'),'linew',2)


ylabel('Normalized max p120ctn int. per cell-pair linescan','Fontsize',15)

[pp120,tblp120,statsp120] = kruskalwallis(AllCellPairs.p120ctnMaxIntNorm,AllCellPairs.tipScoreGroup);
[cp120,mp120,hp120,gnamesp120] = multcompare(statsp120,'Alpha',0.05)






