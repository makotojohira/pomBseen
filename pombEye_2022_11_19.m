% pombEye

% Each column label in the output is preceded by a number and a colon.

% 1: is data from channel 1 (which MUST be a brightfield image)
% 2: is data from channel 2 (which MUST be the main fluorescent image, preferably the experimental channel, often green)
% 3: is optional data from channel 3 (MUST be a fluorescent image, often an internal red control channel)
% 2/2: means channel 2 data masked by nuclei segmented on channel 2
% 2/3: means channel 2 data masked by nuclei segmented on channel 3
% 3/3: means channel 3 data masked by nuclei segmented on channel 3
% 3/2: means channel 3 data masked by nuclei segmented on channel 2


% **MATLAB
% PombEye extensively uses the Matlab image analysis toolbox for many of the
% image manipulation functions, as well as the Bio-Formats toolbox used
% primarily in the dvFileInputs.m function to import the .dv file and split
% them into each channel image and to extract metadata to convert pixel 
% measurements to real units such as micrometers.


% **IMAGES
% PombEye.m works specifically with images taken on a Deltavision 
% fluorescence microscope imaging system and imports .dv image files with
% multiple channels.  However, single-channel bright-field tiff files can 
% also be imported and analyzed for a more limited set of data (cell
% length, width and area).  Data for single-channel images will be reported 
% in pixels, and the user will have to supply their own calibration data 
% and do the conversions.

% Furthermore, this program has been specifically developed and tested on
% S. pombe cells imaged in the bright field mode such that a distinct 
% bright halo defines each cell.  This halo obviously simplifies the 
% segmenting process.

% Image file 2019_09_03_yMO100_27C_1_R3D.dv is an example of a 2-channel
% image that pombEye can process well.

% First step is to execute the import and splitting function from 
% dvFileInputs.m and then exporting the individual .tiff image files.
% dvFileInputs then passes the image and meta data onto subsequent
% functions for processing and analysis.

%%
clc        % Clear the command window
clear      % Clear the workspace
close all  % Close all figure windows

tic        % Start timer
[dvfile,R1,n,R] = dvFileInputs;  % Import Deltavision image file

% dvfile contains both image and metadata information from the Deltavision
% file.  Used bfopen function to extract data, image data assigned to R,
% where image data is then split out into R1. R will be used in
% NuclearCellFilter function, R1 will be used by ImportImage function.

dvFiletime = toc  % Output the time to complete the first function


%%
% Import images, sharpen contrast for BF, convert to grayscale, and return
% 2 or 3 images in a given .dv file.  This ImportImage function contains an
% If/Else operation for single-channel images, in which case the
% ImportImage function runs the full processing and analysis of the single 
% image.

[R1a, R2a, R3a, FN1a, FN2a, FN3a] = ImportImage(dvfile,R1,n);
importImagetime = toc

%%
% Import parameters for various functions throughout pombEye after testing 
% for single-channel image and display message and stop pombEye since all 
% processing and analysis occurred in the ImportImage function.  The
% remainder of pombEye will operate only for multi-channel images.

if n == 1;
    disp('Only one channel data in this image file');
else % everything after this is part of the else statement...
%%
% Import parameters for various functions throughout pombEye after testing 
% for single-channel image.

[ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix, ClearBorderConn, ClearBorderMaxPix, SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax, NucThreshBinNeighbrhd, NucThreshBinMaxPix, NucSegmentConn, NucSegmentMaskMax, NucCellFilterConn, NucCellFilterNCRatioMin, ConvexFilterSlope, ConvexFilterIntercept] = Parameters;


%%
% Apply Otsu's thresholding method, binarize, and return resulting BW
% image for BF channel.
[BW1] = ThreshBinarize(R1a, FN1a, ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix);
threshbintime = toc

%%
% Filter out small regions of noise, invert so cell interiors are white,
% and return clean inverted image.
[BW1b] = InverseBW(BW1, FN1a, InverseBWMaxPix);
invBWtime = toc

%%
% Repeat segment, apply label matrix color code, number cells, and return 
% image.
[CC,Area,BW1b] = SegmentNum(BW1b,FN1a,SegmentNumConn,SegmentNumAreaMin,SegmentNumAreaMax);
segnumtime = toc

%%
% Delete cells contacting image border, repeat filtering, and return image.
[BW1c] = ClearBorder(BW1b,FN1a,ClearBorderConn,ClearBorderMaxPix);
clearbordertime = toc

%%
% Repeat segment, apply label matrix color code, number cells, and return 
% image.
[CC,Area,BW1d] = SegmentNum(BW1c,FN1a,SegmentNumConn,SegmentNumAreaMin,SegmentNumAreaMax);
segnumtime = toc

%%
% Use convex polygon and aspect ratio to filter out non-cell regions that 
% display significant concavity, and return image.
[CCstats,BW1e] = ConvexFilter(CC,BW1d,FN1a,ConvexFilterSlope,ConvexFilterIntercept);
convfiltertime = toc

%%
% Repeat segment, apply label matrix color code, number cells, and return 
% image. 
%BW1d = BW1e;
[CC,Area,BW1f] = SegmentNum(BW1e,FN1a,SegmentNumConn,SegmentNumAreaMin,SegmentNumAreaMax);
segnumtime3 = toc

%%
% Assuming the filtering was not perfect, manually select regions to delete
% and return image.
[BW1f,CC] = SelectRegionDel(BW1f,CC,FN1a);
seldeltime = toc


%%
% End of bright field processing.  Now go on to importing and processing
% fluorescent channel(s).

%%
% First use segmented cells as masks to isolate cellular FITC signal and
% set background to zero - thus eliminating all extraneous signals. Then
% individually segment nuclei from threshold calculated for each cell.  

% NOTES - this works well for bright nuclei - but not as well with
% low-intensity nuclei (small cells).
    
[BW2,THD] = NuclearThreshBinarize(CC,BW1f,R2a,FN2a,NucThreshBinNeighbrhd,NucThreshBinMaxPix);

nucthreshtime = toc
%%
% Find the background intensity for a background subtraction normalization.
% Usually the intensity with the highest count is the background, so run a
% histogram function, and find the maximum count. The bin containing that
% max count should correspond to the background intensity.

% The function imhist is a specific MATLAB function which returns a
% histogram for images.  
[count2a, binloc2a] = imhist(R2a,65536);

% The function max returns the maximum count from that histogram, and the 
% index (which number) of that max count.
[maxCount2a,maxI2a] = max(imhist(R2a,65536));

% So to find the bin corresponding to the max count, I go to the maxI'th
% value in the variable binloc.
maxBin2a = binloc2a(maxI2a);

% The value in maxBin2a is the intensity of the background in image R2a.

%%
% Need to add MaxBin1 to the output spreadsheet
% NOTE - check that the returned maxBin1 is indeed similar to the
% background intensity (e.g., checked with FIJI or similar).

%%
% After binarizing, segment nuclei similar to how cells were segmented -
% but use different size filtering parameters

[CC2, BW2a] = NuclearSegment(BW2,FN2a,NucSegmentConn,NucSegmentMaskMax);
nucsegtime = toc

%%

% Superimpose cell and nuclear segments, count nuclei per cell, and begin
% to analyze cells and nuclei.  Extract calibration data from .dv metadata
% and convert data from pixels to microns.

[CC,cellNucCount] = NuclearCellFilter(BW1f,BW2a,CC,CC2,FN2a,R2a,R,NucCellFilterConn,NucCellFilterNCRatioMin);
nucfiltertime = toc

% Write data table to spreadsheet and save.

varnames1 = {'1: Cell Index','2: Num Nuclei', '1: Cell Length', '2: Nuc Area', '1: Cell Area', '2/2: Mean Nuc Int', '1: Cell Width'};
CellNucTable = array2table(cellNucCount, 'VariableNames',varnames1);
tablename2a = [FN2a(1:end),'_Data.csv'];
writetable(CellNucTable,tablename2a);


%%
% Get whole cell fluorescence and background intensity in ch2:

NucIntensity2 = regionprops(CC,R2a,'MeanIntensity'); % get mean intensity of image in R2a, within masked area in CC

numCells = length(CC.PixelIdxList);
WholeCell2 = [1,numCells];
WholeCell2 = [NucIntensity2.MeanIntensity];
WholeCell2 = WholeCell2';

name8 = '2: Whole Cell Int';
CellNucTable.(name8) = WholeCell2;

maxBinCol2a = [1,numCells];
maxBinCol2a(1:numCells) = maxBin2a;
maxBinCol2a = maxBinCol2a'; 

name9 = '2: Background Int';
CellNucTable.(name9) = maxBinCol2a;

writetable(CellNucTable,tablename2a);

%%
% Subtract background from nuclear and WC fluorescence.

name10 = '2/2: BkgdSubt Nuc Int';
CellNucTable.(name10) = CellNucTable{:,6}-CellNucTable{:,9};


name11 = '2: BkgdSubt WC Int';
CellNucTable.(name11) = CellNucTable{:,8}-CellNucTable{:,9};

writetable(CellNucTable,tablename2a);

%%
% End of 2nd channel fluorescent image processing.  This next section is
% optional, only if a 3rd channel fluorescent image is taken.

%%
% Test if there is a 3rd channel, if not then display message, but if so
% then proceed to 3rd channel image processing and analysis.

if isempty(R3a)
    disp('No third channel data in this image file');
else   % everything after this is part of the else statement...


%%
% Import the TRITC image, convert to grayscale, and return image.
% Try to make this only occur if there is three channels of data in the .dv
% file.

% First use segmented cells as masks to isolate cellular TRITC signal and
% set background to zero - thus eliminating all extraneous signals. Then
% individually segment nuclei from threshold calculated for each cell.  

[BW3,THD] = NuclearThreshBinarize(CC,BW1f,R3a,FN3a,NucThreshBinNeighbrhd,NucThreshBinMaxPix);
nucthreshtime = toc

%%
% Find the background intensity so user might do a background subtraction.
% Usually the intensity with the highest count is the background, so run a
% histogram function, and find the maximum count. The bin containing that
% max count should correspond to the background intensity.

% The function imhist is a specific MATLAB function which returns a
% histogram for images.  
[count3a, binloc3a] = imhist(R3a,65536);

% The function max returns the maximum count from that histogram, and the 
% index (which number) of that max count.
[maxCount3a,maxI3a] = max(imhist(R3a,65536));

% So to find the bin corresponding to the max count, I go to the maxI'th
% value in the variable binloc.
maxBin3a = binloc3a(maxI3a);

% The value in maxBin2a is the intensity of the background in image R2a.

%%
% After binarizing, segment nuclei similar to how cells were segmented -
% but use different size filtering parameters

[CC3, BW3a] = NuclearSegment(BW3,FN3a,NucSegmentConn,NucSegmentMaskMax);
nucsegtime = toc

%%
% Next is NuclearCellFilter as before for 2nd channel, but ch3 fluorescent
% data masked by segmented ch2 nuclei...
[CC,cellNucCount3] = NuclearCellFilter(BW1f,BW2a,CC,CC2,FN3a,R3a,R,NucCellFilterConn,NucCellFilterNCRatioMin);
nucfiltertime = toc

name12 = '3/2: Mean Nuc Int';
CellNucTable.(name12) = cellNucCount3(:,6);
%tablename3b = [FN3a(1:end),'_Data_1.csv'];
%writetable(CellNucTable,tablename3b);
writetable(CellNucTable,tablename2a);


%%
% Now get channel 3 data using nuclei segmented from ch3 (previous 
% NuclearCellFilter just above analyzed ch3 fluorescence using nuclei 
% segmented from ch2):

[CC,cellNucCount3a] = NuclearCellFilter(BW1f,BW3a,CC,CC3,FN3a,R3a,R,NucCellFilterConn,NucCellFilterNCRatioMin);
nucfiltertime = toc

name13  = '3: Nuc Area';
CellNucTable.(name13) = cellNucCount3a(:,4);

name14 = '3/3: Mean Nuc Int';
CellNucTable.(name14) = cellNucCount3a(:,6);

writetable(CellNucTable,tablename2a);


%%
% Now get channel 2 data using nuclei segmented from ch3:

[CC,cellNucCount3a] = NuclearCellFilter(BW1f,BW3a,CC,CC3,FN3a,R2a,R,NucCellFilterConn,NucCellFilterNCRatioMin);
nucfiltertime = toc

name15 = '2/3: Mean Nuc Int';
CellNucTable.(name15) = cellNucCount3a(:,6);

writetable(CellNucTable,tablename2a);



%%
% Get whole cell and background intensity fluorescence in ch3:

NucIntensity3 = regionprops(CC,R3a,'MeanIntensity'); % get mean intensity of image in R3a, within masked area in CC

WholeCell3 = [1,numCells];
WholeCell3 = [NucIntensity3.MeanIntensity];
WholeCell3 = WholeCell3';

name16 = '3: Whole Cell Int';
CellNucTable.(name16) = WholeCell3;

maxBinCol3a = [1,numCells];
maxBinCol3a(1:numCells) = maxBin3a;
maxBinCol3a = maxBinCol3a'; 

name17 = '3: Background Int';
CellNucTable.(name17) = maxBinCol3a;

writetable(CellNucTable,tablename2a);

%%
% Subtract background from nuclear and WC fluorescence.

name18 = '3/3: BkgdSubt Nuc Int';
CellNucTable.(name18) = CellNucTable{:,14}-CellNucTable{:,17};

name19 = '3: BkgdSubt WC Int';
CellNucTable.(name19) = CellNucTable{:,16}-CellNucTable{:,17};

name20 = '3/2: BkgdSubt Nuc Int';
CellNucTable.(name20) = CellNucTable{:,12}-CellNucTable{:,17};

name21 = '2/3: BkgdSubt Nuc Int';
CellNucTable.(name21) = CellNucTable{:,15}-CellNucTable{:,17};



writetable(CellNucTable,tablename2a);

%%
% Add a few more parameters to add to the data table

% number of nuclei/cell in channel 3
name22 = '3: Num Nuclei';
CellNucTable.(name22) = cellNucCount3a(:,2);

% width of cells in channel 3
% name15 = '3: Cell Width';
% CellNucTable.(name15) = cellNucCount3a(:,7);

writetable(CellNucTable,tablename2a);

end
end



