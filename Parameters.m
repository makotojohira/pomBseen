% Set up a parameters page so all hard-coded parameters for
% various functions and calculations are compiled and called from one
% place.


function [ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix, ClearBorderConn, ClearBorderMaxPix, SegmentNumConn, SegmentNumAreaMin, SegmentNumAreaMax, NucThreshBinNeighbrhd, NucThreshBinMaxPix, NucSegmentConn, NucSegmentMaskMax, NucCellFilterConn, NucCellFilterNCRatioMin, ConvexFilterSlope, ConvexFilterIntercept] = Parameters
%ConvexFilterSlope = 12.8571
%ConvexFilterIntercept = 12.5


ThreshBinSensitivity = 0.55;
% In function ThreshBinarize, std function adaptthresh has a sensitivity 
% parameter I have set at 0.55

ThreshBinNeighbrhd = [15 15];
% In function ThreshBinarize, std function adaptthresh has a 
% neighborhoodsize parameter I have set at [15 15]

InverseBWMaxPix = 600;
% In function InverseBW, std function bwareaopen has a maximum pixels 
% parameter I have set at 600; this function is called before and after the 
% BW inversion

ClearBorderConn = 8;
% In function ClearBorder, std function imclearborder has a maximum pixels 
% parameter I have set at 8

ClearBorderMaxPix = 1200;
% In function ClearBorder, std function bwareaopen has a connectivity 
% parameter I have set at 1200

SegmentNumConn = 4;
% In function SegmentNum, std function bwconncomp has a connectivity 
% parameter I have set at 4

SegmentNumAreaMin = 500;
SegmentNumAreaMax = 100000;
% In function SegmentNum, a mask I set up has a area parameter 
% <500 | >100000

%ConvexFilter Eqn: y = 12.8571x - 12.5
% In function ConvexFilter, a mask I set up has an equation to retail cells 
% and exclude artifacts
ConvexFilterSlope = 12.8571;
ConvexFilterIntercept = 12.5;

NucThreshBinNeighbrhd = ones(5,5);
% In function NuclearThreshBinarize, std function imopen has a neighborhood 
% parameter I have set at ones(5,5)

NucThreshBinMaxPix = 200;
% In function NuclearThreshBinarize, std function Bbwareaopen has a maximum 
% pixels parameter I have set at 200

NucSegmentConn = 4;
% In function NuclearSegment, std function bwconncomp has a connectivity 
% parameter I have set at 4

NucSegmentMaskMax = 155;
% In function NuclearSegment, a mask I set up has a parameter NucArea I 
% have set at <155

NucCellFilterConn = 4;
% In function NuclearCellFilter, std function bwconncomp has a connectivity 
% parameter I have set at 4

NucCellFilterNCRatioMin = 0.2;
% In function NuclearCellFilter, a mask I have set up has a parameter 
% NucArea(end)/CellArea(end) > 0.2




end
