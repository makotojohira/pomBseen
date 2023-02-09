% pombEyeLength

% 2022 10 25 - modify to be called within pombEye as a function in an IF
% statement.

% Revision 2021_11_30 was modified for Expt 101J on 2022 08 02 in order to
% analyze BF images only and return cell size.  This has been tested on
% images from the Rhind Lab's Zeiss fluorescence scope.

% NOTE - the length conversion is based on Samir's 40x calibration of 7.5
% pixels per micron (and confirmed by MO).  If you take a 100X image,
% change the conversion from 7.5 to 18.7 in line 169 of this code.

% **IMAGES

% This program has been specifically developed and tested on
% S. pombe cells imaged in the bright field (BF) mode such that a distinct 
% bright halo defined each cell.  This halo obviously simplifies the 
% segmenting process.

% Image file 2022 07 15 Expt 101J N1-0001.tif is an example of a BF
% image that this pombEye rev can process well.

% **MATLAB
% PombeY extensively uses the Matlab image analysis toolbox for many of the
% image manipulation functions, as well as the Bio-Formats toolbox used
% primarily in the dvFileInputs.m function to import the .dv file and split
% them into each channel image and to extract metadata to convert pixel 
% measurements to real units such as micrometers.

% First step is to execute the import and splitting function from 
% dvFileInputs.m and then exporting the individual .tiff image files.

%%

function [] = pombEyeLength(dvfile,R1,n,R1a,FN1a)

% Clear the command window
%clc        

% Clear the workspace
%clear      

% Close all figure windows
%close all  

% Start timer
%tic        

% Import Deltavision image file
% [dvfile,R1,R] = dvFileInputs;  Already done in pombEye

% dvfile contains both image and metadata information from the Deltavision
% file.  Used bfopen function to extract data, image data assigned to R,
% where image data is then split out into R1. R will be used in
% NuclearCellFilter function, R1 will be used by ImportImage function.

% Output the time to complete the first function
%dvFiletime = toc  

%%

% Import images, sharpen contrast for BF, convert to grayscale, and return
% 2 or 3 images in a given .dv file.
%[R1a, R2a, R3a FN1a FN2a FN3a] = ImportImage(dvfile,R1,n);
%[R1a, FN1a] = ImportImage(dvfile,R1,n);

% Output the time to complete the second function
%pELImportImagetime = toc


%clearvars -except R1a FN1a
%%
    %%
    % start copy-paste from pombEyeLength

% Import parameters for various functions throughout pombEye.

[ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix, ClearBorderConn, ClearBorderMaxPix, SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax, NucThreshBinNeighbrhd, NucThreshBinMaxPix, NucSegmentConn, NucSegmentMaskMax, NucCellFilterConn, NucCellFilterNCRatioMin, ConvexFilterSlope, ConvexFilterIntercept] = Parameters;


%%

% Apply Otsu's thresholding method, binarize, and return resulting BW
% image for BF channel.
[BW1] = ThreshBinarize(R1a, FN1a);

% Output the time to complete the third function
pELthreshbintime = toc

%%

% Filter out small regions of noise, invert so cell interiors are white,
% and return clean inverted image.
[BW1b] = InverseBW(BW1, FN1a);

% Output the time to complete the 4th function
pELinvBWtime = toc

%%

% Delete cells contacting image border, repeat filtering, and return image.
[BW1c] = ClearBorder(BW1b,FN1a);

% Output the time to complete the 5th function
pELclearbordertime = toc

%%

% Repeat segment, apply label matrix color code, number cells, and return 
% image.
[CC,Area,BW1d] = SegmentNum(BW1c,FN1a);

% Output the time to complete the 6th function
pELsegnumtime = toc

%%

% Use convex polygon and aspect ratio to filter out non-cell regions that 
% display significant concavity, and return image.
[CCstats,BW1e] = ConvexFilter(CC,BW1d,FN1a);

% Output the time to complete the 7th function
pELconvfiltertime = toc

%%

% Repeat segment, apply label matrix color code, number cells, and return 
% image. 
%BW1d = BW1e;
[CC,Area,BW1f] = SegmentNum(BW1e,FN1a);

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

Length = regionprops(CC, 'MajorAxisLength');

cellNucCount = [numCells,1];

cellNucCount(1:numCells,1) = 1:numCells;

cellNucCount(1:numCells,2) = [Length.MajorAxisLength];

disp('Length');

% Note that these Length values are in pixels - we want to convert that to
% microns.  Samir did a calibration and I will use his values, used in
% previous calculations for this project.

cellNucCount(1:numCells,3) = [Length.MajorAxisLength]/7.5;

disp('cellNucCount');


%%
% Now I want to report lengths:
% Write the data in cellNucCount to a table and save it
CellNucTable = array2table(cellNucCount, 'VariableNames',{'Cell Index','Cell Length (pixels)','Cell Length (microns)'});
tablename2a = [FN1a(1:end),'_Data.csv'];
writetable(CellNucTable,tablename2a);

disp('writetable');

%%

clearvars 

%%
% end copy-paste from pombEyeLength





