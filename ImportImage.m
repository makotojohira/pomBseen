function [R1a R2a R3a FN1a FN2a FN3a] = ImportImage(dvfile,R1,n)

% Import BF image, convert to grayscale, and return/save image.
% The 1st image must be BF.  The 2nd and 3rd must be fluorescent images.  
% Segment on the BF image, specifically using the halo around each cell.

% Since we are specifying microscopy images in which cells are reliably
% surrounded by a bright halo, we will sharpen the image to emphasize that
% halo which we will use in later steps to segment each cell.  Will not do
% that for the fluorescent images.

R1a = []; 
R2a = [];
R3a = [];
FN1a = [];
FN2a = [];
FN3a = [];

R1a = imsharpen(R1{1,1});
R1a = mat2gray(R1a);
FN1a = [dvfile(1:end-3) '_1'];
figure('Numbertitle', 'off','Name','Function: pELImportImage');
imshow(R1a);
pause(1);
title(FN1a, 'Interpreter', 'none');


%R2a = mat2gray(R1{2,1});


if n == 1;
    disp('Only one channel data in this image file');
    %%
    %%
    % start copy-paste from pombEyeLength
% Import parameters for various functions throughout pombEye.

[ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix, ClearBorderConn, ClearBorderMaxPix, SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax, NucThreshBinNeighbrhd, NucThreshBinMaxPix, NucSegmentConn, NucSegmentMaskMax, NucCellFilterConn, NucCellFilterNCRatioMin, ConvexFilterSlope, ConvexFilterIntercept] = Parameters;


%%

% Apply Otsu's thresholding method, binarize, and return resulting BW
% image for BF channel.
[BW1] = ThreshBinarize(R1a, FN1a, ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix);

% Output the time to complete the third function
pELthreshbintime = toc

%%

% Filter out small regions of noise, invert so cell interiors are white,
% and return clean inverted image.
[BW1b] = InverseBW(BW1, FN1a, InverseBWMaxPix);

% Output the time to complete the 4th function
pELinvBWtime = toc

%%

% Repeat segment, apply label matrix color code, number cells, and return 
% image.
[CC,Area,BW1b] = SegmentNum(BW1b,FN1a,SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax);

% Output the time to complete the 6th function
pELsegnumtime = toc

%%

% Delete cells contacting image border, repeat filtering, and return image.
[BW1c] = ClearBorder(BW1b,FN1a,ClearBorderConn, ClearBorderMaxPix);

% Output the time to complete the 5th function
pELclearbordertime = toc

%%

% Repeat segment, apply label matrix color code, number cells, and return 
% image.
[CC,Area,BW1d] = SegmentNum(BW1c,FN1a,SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax);

% Output the time to complete the 6th function
pELsegnumtime = toc

%%

% Use convex polygon and aspect ratio to filter out non-cell regions that 
% display significant concavity, and return image.
[CCstats,BW1e] = ConvexFilter(CC,BW1d,FN1a,ConvexFilterSlope, ConvexFilterIntercept);

% Output the time to complete the 7th function
pELconvfiltertime = toc

%%

% Repeat segment, apply label matrix color code, number cells, and return 
% image. 
%BW1d = BW1e;
[CC,Area,BW1f] = SegmentNum(BW1e,FN1a,SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax);

% Output the time to complete the 8th function
pELsegnumtime3 = toc

%%
% 2022 08 02
% Now extract and modify code from NuclearCellFilter to find cell lengths:

CC = bwconncomp(BW1f,4);
s = regionprops(CC,'centroid');
centroid = cat(1,s.Centroid);
n = CC.NumObjects;

overlay = imfuse(BW1f, BW1e);
figure('Numbertitle', 'off','Name','Function: NuclearFilter.m Overlay');;
imshow(overlay);
%nucfilter1 = toc(nucfilterstart)

for n=1:n;
    text(centroid(n,1),centroid(n,2),sprintf('%d',n),'HorizontalAlignment','center');
end

drawnow
hold off

disp('Figure overlay');

% Now count nuclear segments in each cell.  NOTE - I don't want to do
% this...

% First, how many cells?  Another way to do this
numCells = length(CC.PixelIdxList);

disp('Numcells');

% Do a for loop, and for each cell see if I can count the number of nuclear
% segments - Comment out
%labeled2 = labelmatrix(CC2);
%labeled2 = double(labeled2);
%cellNucCount = [numCells,2];
%nucfilter2 = toc(nucfilterstart)

% 2022 08 02 - extract length code from for loop:

Length = regionprops(CC, 'MajorAxisLength', 'MinorAxisLength', 'Area');

cellNucCount = [numCells,1];

cellNucCount(1:numCells,1) = 1:numCells;
cellNucCount(1:numCells,2) = [Length.MajorAxisLength];
cellNucCount(1:numCells,3) = [Length.MinorAxisLength];
cellNucCount(1:numCells,4) = [Length.Area];


disp('Length');

% Note that these Length values are in pixels - we want to convert that to
% microns.  Samir did a calibration and I will use his values, used in
% previous calculations for this project.

%cellNucCount(1:numCells,5) = [Length.MajorAxisLength]/7.5;

disp('cellNucCount');


%%
% Now I want to report lengths:
% Write the data in cellNucCount to a table and save it
% CellNucTable = array2table(cellNucCount);
% tablename2a = [FN1a(1:end),'_Data.csv'];
% CellNucTable = table(CellNucTable, 'VariableNames',{'Cell_Index','Cell_Length_pixels','Cell_Length_microns'});
% writetable(CellNucTable,tablename2a);

%CellNucTable = array2table(cellNucCount, 'VariableNames',{'Cell_Index','Cell_Length_pixels','Cell_Length_microns'});

varNames = {'1: Cell Index','1: Cell Length pixels', '1: Cell Width pixels', '1: Cell area pixels^2'};  %,'1: Cell Length microns'};
CellNucTable = array2table(cellNucCount, 'VariableNames',varNames);
tablename2a = [FN1a(1:end),'_Data.csv'];
writetable(CellNucTable,tablename2a);

disp('writetable');


%%
% end copy-paste from pombEyeLength
%%
%%

else   % everything after this is part of the else statement...

R2a = R1{2,1};
FN2a = [dvfile(1:end-3) '_2'];
figure('Numbertitle', 'off','Name','Function: ImportImage');
imshow(mat2gray(R2a));
title(FN2a, 'Interpreter', 'none');


FN3a = [];
R3a = [];
if n == 3
    R3a = R1{3,1};
    FN3a = [dvfile(1:end-3) '_3'];
    figure('Numbertitle', 'off','Name','Function: ImportImage');
    imshow(mat2gray(R3a));
    title(FN3a, 'Interpreter', 'none');
end

end

clearvars -except R1 R1a R2a R3a FN1a FN2a FN3a
