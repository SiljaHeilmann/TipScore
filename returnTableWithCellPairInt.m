function TcellPairs = returnTableWithCellPairInt(NUCLEI_BW,IM,approxPixelWidthOfCellNuclei,plot_yes_no)

TcellPairs = table;% initialize table for storing data on cells

L = bwlabel(NUCLEI_BW,4);

r = regionprops(L,'centroid');

centroids = cat(1,r.Centroid);

D = squareform(pdist(centroids));% a pairwise distance matrix. D(i,j) holds the distance between the centroid in nucleus i and j

Dbw = D<approxPixelWidthOfCellNuclei.*1.5;% locate neiborpairs (cell pairs whos centroids are less than approxPixelWidthOfCellNuclei.*1.5 appart)
Dbw(D==0)=0; % Diagonal (self pairs) should be zero
Dbw = triu(Dbw);% number of ones in Dbw is the number of cell neigbor pair that we will look at (upper triangle of matrix is enough since its symmetric)

TcellPairs.cellPairIDnumber = [1:sum(Dbw(:))]';% assign each cell pair a number

% initialise variable columns
TcellPairs.IDcell1  = NaN(size([1:sum(Dbw(:))]'));
TcellPairs.IDcell2  = NaN(size([1:sum(Dbw(:))]'));
TcellPairs.coorx    = NaN(size([1:sum(Dbw(:))]'));
TcellPairs.coory    = NaN(size([1:sum(Dbw(:))]'));

TcellPairs.lineScan       = cell(size([1:sum(Dbw(:))]'));
TcellPairs.lineScanLength = NaN(size([1:sum(Dbw(:))]'));
TcellPairs.lineScanMax    = NaN(size([1:sum(Dbw(:))]'));
TcellPairs.lineScanMed    = NaN(size([1:sum(Dbw(:))]'));


IMmemb = IM;

% gauss filt membrane channel image (to reduce noise on a scale much smaller than a nucleus)
IMmemb = imgaussfilt(IMmemb,[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(0.75/50));

if plot_yes_no == 1
    
    figure('units','normalized','outerposition',[0 0.2 1 0.6]);
    imshow(double(IMmemb)./prctile(double(IMmemb(:)),95),[]); hold on
    plot(centroids(:,1),centroids(:,2),'yo')
end
linescanMatrix = NaN(sum(Dbw(:)),ceil(approxPixelWidthOfCellNuclei.*1.5)+1);% number of cell neigbor pairs x maximun length of line scan
% tipscore =  NaN(sum(Dbw(:)),1);% number of cell neigbor pairs

cc = 1; % cell pair counter
for ii=1:max(L(:))% iterate through all cell
    indexCurrentNeigbors = find(Dbw(ii,:));% find index of current cells neigbors
    for jj=1:length(indexCurrentNeigbors)% iterate through each neigbor
            
        % get pixel location of cell pair to extract tip score for current
        % cell pair
        middlepoint_x = round(mean([centroids(ii,1) centroids(indexCurrentNeigbors(jj),1)]));
        middlepoint_y = round(mean([centroids(ii,2) centroids(indexCurrentNeigbors(jj),2)]));
        
        if plot_yes_no ==1
            plot([centroids(ii,1) centroids(indexCurrentNeigbors(jj),1)],[centroids(ii,2) centroids(indexCurrentNeigbors(jj),2)],'-w','LineWidth',1)
            %plot(middlepoint_x,middlepoint_y,'.r','MarkerSize',10); hold on % plot point where tipscore is extracted
        end        
        % store tipscore of current cell pairs
        %   tipscore(cc) = TSC(middlepoint_y,middlepoint_x); % note this may not be the exact location of the max pixel along the line scan but it will be very close and is good enough...
        
        % Extract line scan between cell centroids
        val =  improfile(IMmemb,[centroids(ii,1) centroids(indexCurrentNeigbors(jj),1)],[centroids(ii,2) centroids(indexCurrentNeigbors(jj),2)]);
        linescanMatrix(cc,1:length(val)) = val;% store linescan in matrix
        
        TcellPairs.coorx(TcellPairs.cellPairIDnumber==cc) = middlepoint_x;
        TcellPairs.coory(TcellPairs.cellPairIDnumber==cc) = middlepoint_y;
        TcellPairs.IDcell1(TcellPairs.cellPairIDnumber==cc) = ii;
        TcellPairs.IDcell2(TcellPairs.cellPairIDnumber==cc) = indexCurrentNeigbors(jj);
        
        TcellPairs.lineScan(TcellPairs.cellPairIDnumber==cc) = {val};
        TcellPairs.lineScanLength(TcellPairs.cellPairIDnumber==cc) = length(val);
        TcellPairs.lineScanMax(TcellPairs.cellPairIDnumber==cc) = max(smoothdata(val,1,'movmean',3));
        TcellPairs.lineScanMed(TcellPairs.cellPairIDnumber==cc) = median(val);
        
        cc = cc + 1;% go to next cell pair
        
    end
end

end

