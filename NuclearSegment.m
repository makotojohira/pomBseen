
% After binarizing now segment nuclei similar to how cells were segmented -
% but use different size filtering parameters

% Segment, apply label matrix color code, and return image for nuclei.

function [CC2, BW2a] = NuclearSegment(BW2,FN2a,NucSegmentConn,NucSegmentMaskMax)

% Modified stats.Solidity from <0.90 to <0.25 and that made
% a huge difference in picking up larger cells - but this picked up
% branched cells so I set it at <0.5 and this seemed to find a happy
% balance.  I also changed the mask to Area <150 and >a ratio (see below).

CC2      = bwconncomp(BW2,NucSegmentConn); 
% NucSegmentConn = 4
stats   = regionprops(CC2,{'Area' 'Solidity'});
NucArea    = [stats.Area];
mask    = NucArea < NucSegmentMaskMax; 
% NucSegmentMaskMax = 155 defined in the Parameters function

CC2.PixelIdxList(mask) = [];
CC2.NumObjects   = length(CC2.PixelIdxList);
NucArea(mask)      = [];
BW2a      = false(size(BW2));
BW2a(vertcat(CC2.PixelIdxList{:})) = true;

% Now to visualize the segmenting step above:
s = regionprops(CC2,'centroid');
centroid = cat(1,s.Centroid);
n = CC2.NumObjects;

labeled = labelmatrix(CC2);
RGB_label2 = label2rgb(labeled,'spring','c','shuffle');

figure('Numbertitle', 'off','Name','Function: NuclearSegment.m');
imshow(RGB_label2);
pause(1);
hold on;
%plot(centroid(:,1),centroid(:,2),'b*') % Confirm location of centroids
for n=1:n;
    text(centroid(n,1),centroid(n,2),sprintf('%d',n),'HorizontalAlignment','center');
end
hold off;
title(FN2a, 'Interpreter', 'none');

clearvars -except CC2 BW2a

% Use label2rgb to choose the colormap, the background color, and how 
% objects in the label matrix map to colors in the colormap. In the 
% pseudocolor image, the label identifying each object in the label matrix 
% maps to a different color in an associated colormap matrix.
