
% Filtering for nuclear/cell area - used 0.2 as cuttoff - this cutoff
% may be modified.

% Then filter out cells with either < or > than 1 nuclei per cell.  Compare 
% the segmented cell matrix, and the segmented nuclei matrix, and set the 
% cells to zero if it overlaps with either 0 or >1 nuclei...

function [CC,cellNucCount] = NuclearCellFilter(BW1f,BW2a,CC,CC2,FN2a,R2a,R,NucCellFilterConn,NucCellFilterNCRatioMin);

nucfilterstart = tic

% First let's superimpose the images of previous fully filtered cells and
% segmented nuclei

CC = bwconncomp(BW1f,NucCellFilterConn);  
% NucCellFilterConn = 4 defined in the Parameters function

s = regionprops(CC,'centroid');
centroid = cat(1,s.Centroid);
n = CC.NumObjects;

overlay = imfuse(BW1f, BW2a);
figure('Numbertitle', 'off','Name','Function: NuclearFilter.m Overlay');;
imshow(overlay);
nucfilter1 = toc(nucfilterstart)
pause(1);
for n=1:n;
    text(centroid(n,1),centroid(n,2),sprintf('%d',n),'HorizontalAlignment','center');
end

drawnow
hold off


% Now count nuclear segments in each cell.

% First, how many cells?  Another way to do this

numCells = length(CC.PixelIdxList);

% Do a for loop, and for each cell see if I can count the number of nuclear
% segments
labeled2 = labelmatrix(CC2);
labeled2 = double(labeled2);
cellNucCount = [numCells,2];
nucfilter2 = toc(nucfilterstart)

%blank = 0;

for n = 1:numCells
    cell = false(size(BW1f));   % each time blank the image area
    cell(CC.PixelIdxList{n}) = true;   % each time display a new cell by index
    %cell = double(cell);   % change type so we can multiply
    %cellNuc = cell.* labeled2;   % overlap each cell with the associated nuclei that fall within its area
    cellNuc = labeled2;  % start filtering by assigning labeled2 to new working var cellNuc
    cellNuc(~cell) = 0;  % now remove any nuclear signal not associated with indexed cells
    %figure;
    %imshow(cellNuc);
    count = unique(cellNuc);   % find the unique values which will always include 0 for background, and a value for each nuclei
    nNuc = length(count) - 1;
    cellNucCount(n,1) = n;   % create a table and assign each cell index
    cellNucCount(n,2) = nNuc;   % and associated with each cell index is the number of nuclei
    Nindex = cellNucCount(n,2) < 1 | cellNucCount(n,2) > 1;  % logical index true if < or > 1, false if 1
    % Add length of individual cell to table cellNucCount
    Length = regionprops(cell, 'MajorAxisLength');
    cellNucCount(n,3) = [Length.MajorAxisLength];
    % Now add nuclear intensity to same table if only a single nucleus
    % otherwise enter "N/A".  Also, only if area ratio <0.2
    NucArea = [];  % blank this variable each iteration
    if nNuc == 0
        NucArea = 1;
    else
        NucArea = regionprops(cellNuc, 'Area');  % Get the mean nuclear area from the filtered set of nuclei
        NucArea = [NucArea.Area];
    end;
    CellArea = []; % blank this variable each iteration
    CellArea = regionprops(cell, 'Area');  % Get mean cell area
    CellArea = [CellArea.Area];
    cellNucCount(n,4) = NucArea(end);
    cellNucCount(n,5) = CellArea;
    NucIntensity = []; % blank this variable each iteration
    NCRatioIndex = NucArea(end) / CellArea > NucCellFilterNCRatioMin;  % Index of cells with excessive nuclear size
            % NucCellFilterNCRatioMin = 0.2 defined in the Parameters function
    if NCRatioIndex == 1 | nNuc == 0 %| nNuc ~= 1
        %cellNucCount(n,6) = 0;
        NucIntensity = 0;
    else
        NucIntensity = regionprops(cellNuc,R2a,'MeanIntensity'); % get mean intensity of image in R2a, within masked area in cellNuc
        NucIntensity = [NucIntensity.MeanIntensity];
    end
    cellNucCount(n,6) = NucIntensity(end);
    %    cellNucCount(n,6) = [NucIntensity(end).MeanIntensity];
    Width = regionprops(cell, 'MinorAxisLength');
    cellNucCount(n,7) = [Width.MinorAxisLength];
    clear NucArea; % blank this variable each iteration
end
nucfilter3 = toc(nucfilterstart)


%%
% Here is code copied from NuclearCellData (2021-04-02):

omeMeta = R{1,4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value;
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value;

% need to convert a java data type to a Matlab type
voxelSizeX = double(voxelSizeX);

cellNucCount(:,3) = cellNucCount(:,3) * voxelSizeX;
cellNucCount(:,4) = cellNucCount(:,4) * (voxelSizeX)^2;
cellNucCount(:,5) = cellNucCount(:,5) * (voxelSizeX)^2;
cellNucCount(:,7) = cellNucCount(:,7) * voxelSizeX;

%%

clearvars -except CC cellNucCount

% We also need to delete the double nuclei for those cells that were
% deleted in the previous step - either that or make sure that only the
% single nuclei associated with the fully filtered cells are measured and
% plotted









