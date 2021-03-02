function Tcells = returnTableWithCellInt(NUCLEI_BW,IM,approxPixelWidthOfCellNuclei)
% returnTableWithCellInt(): This function takes a 2D nuclei BW mask and an
% 2D image and returns a table that contains the image intensity in each
% nucleus centroid

% INPUT:
% NUCLEI_BW: 2D BW mask (1's where there are nuclei and 0s elsewhere) of type, logical, double, uint8 or uint16
% IM: 2D Intensity image of type double, uint8 or uint16. Could eg. be a
% immunohistochemical staining of a nuclear marker.
% approxPixelWidthOfCellNuclei: Length scale/diameter of an average
% nucleus. Used for setting the width of the gaussian filtering of the
% image. (We filter/blur the image to ensure that the extracted pixel
% intensities are representative of the local intensities and not extreme outlyers. If you dont want any gaussian filtering then set 'approxPixelWidthOfCellNuclei=1').
%
% OUTPUT:
% Tcells: table object with 4 columns and one row for each cell/nucleus.
% Col 1 have index number/ ID number in labelled image
% Col 2 has x coordinate of nucleus centroid
% Col 3 has y coordinate of nucleus centroid
% Col 4 has intensity value in gaussian filtered image IM at nucleus
% centroid location

Tcells = table;% initialize table for storing data on cells

L = bwlabel(NUCLEI_BW,4);

Tcells.CellIDnumber = [1:max(L(:))]';% 1st column gets ID numbers

r = regionprops(L,'centroid');

centroids = cat(1,r.Centroid);

Tcells.coorx = ceil(centroids(:,1));% 2st column gets xcoor
Tcells.coory = ceil(centroids(:,2));% 3st column gets ycoor

IMorg = IM;

IM(NUCLEI_BW==0) = NaN;% set values outside nuclei mask to NaN in order to avoid that signal from outside values in cytoplasm after gauss filtering
IM = imgaussfilt(double(IM),[approxPixelWidthOfCellNuclei approxPixelWidthOfCellNuclei].*(3/50));

ind = sub2ind(size(IM),Tcells.coory,Tcells.coorx);% return liniar index for centroid pixel locations

Tcells.IMintensity = IM(ind);

% figure
% %imshowpair(IM./max(IM(:),NUCLEI_BW),[0 0.1]); hold on
% %imshowpair(5.*uint8(255.*(IM./max(IM(:)))),bwperim(NUCLEI_BW)); hold on
% imshowpair(5.*uint8(255.*(IM./max(IM(:)))),5.*uint8(255.*(IMorg./max(IMorg(:))))); hold on
% 
% plot(Tcells.coorx,Tcells.coory,'yo'); title('Intensity image and nuclei locations')

end

