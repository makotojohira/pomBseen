
% Tried tophat filtering and single-threshold methods - but went back to
% the for loop.  I did implement post-treshold morphological operators like
% imfill and imopen.

% First use segmented cells as masks to isolate cellular FITC signal and
% set background to zero - thus eliminating all extraneous signals. Then
% use multithresholding to set cell background to zero and thus isolate
% nuclear signal above background.  But found a way to individually segment
% nuclei from threshold calculated for each cell.  Probably not most
% efficient way to do this but used a for loop to do this cell by cell.

% NOTES - this works well for bright nuclei - but not so well with
% low-intensity nuclei (small cells).

function [BW2,THD] = NuclearThreshBinarize(CC,BW1e,R2a,FN2a,NucThreshBinNeighbrhd,NucThreshBinMaxPix)

nucthreshstart = tic

figure;
imshow(mat2gray(R2a));

nuc1time = toc(nucthreshstart)

numCells = length(CC.PixelIdxList);

R2a = mat2gray(R2a);  % Convert R2a to grayscale from 0 to 1 range

% Create a matrix of threshold values for each cell
BW1 = zeros(size(R2a));
BW2 = zeros(size(R2a));
centroidlist1(1:numCells,2) = zeros;

for n = 1:numCells;
    cell = false(size(BW1e));  % each time blank the image area
    cell(CC.PixelIdxList{n}) = true;  % each time display a new cell by index
    R2maskN = R2a;  % assign FITC image to a working variable
    R2maskN(~cell) = 0;  % keep only the FITC data within the indexed cell
    R2maskN1 = R2maskN(R2maskN > 0); % Only consider individual masked cell area for thresholding
    THD = graythresh(R2maskN1);  % Use Otsu's method to threshold the nuclei within the cell perimeter
    BW2 = imbinarize(R2maskN, THD);  % Now apply threshold to the FITC image of the cell - this should identify the nuclei within the cell
    BW2 = imfill(BW2, 'holes');
    BW2 = imopen(BW2, NucThreshBinNeighbrhd);  % NucThreshBinNeighbrhd = one(5,5)
    BW2 = bwareaopen(BW2, NucThreshBinMaxPix); 
    % NucThreshBinMaxPix = 200 defined in the Parameters function
    
    BW1 = BW1 + BW2;  % Build up the BW image of all the nuclei one by one
    % s2 = regionprops(BW2,'centroid');  % Obtain region properties (centroid) for the nuclear BW thresholded image
    % centroidlist1(n,:) = s2.Centroid;  % Build list of centroid X,Y locations
end
nuc2time = toc(nucthreshstart)

BW2 = BW1;

% Now I want to see what this looks like...
%figure('Numbertitle', 'off','Name','Function: NuclearThreshBinarize.m');
%imshow(BW2);
%title(FN2a, 'Interpreter', 'none');
%nuc3time = toc(nucthreshstart)

% BW2 = imfill(BW2, 'holes');
% BW2 = imopen(BW2, ones(5,5));
% BW2 = bwareaopen(BW2, 200);
figure('Numbertitle', 'off','Name','Function: NuclearThreshBinarize.m - cleanup');
imshow(BW2);
pause(1);
hold on;
%for n=1:numCells;
%    text(centroidlist1(n,1),centroidlist1(n,2),sprintf('%d',n),'HorizontalAlignment','center');
%end
title(FN2a, 'Interpreter', 'none');
%drawnow;
hold off;
nuc4time = toc(nucthreshstart)


clearvars -except BW2 THD

