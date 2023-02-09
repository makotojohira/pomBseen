
% Filter out small regions of noise, invert so cell interiors are white,
% and return clean inverted image.
function [BW1b] = InverseBW(BW1, FN1a, InverseBWMaxPix)

% remove small areas less than 600 pixels in area in the cells

bwstart = tic
BW1a = bwareaopen(BW1,InverseBWMaxPix);
% InverseBWMaxPix = 600 defined in the Parameters function

bw1time = toc(bwstart)

% now reverse the image so cells are white

BW1b = imcomplement(BW1a);
bw2time = toc(bwstart)

figure('Numbertitle', 'off','Name','Function: InverseBW.m');
imshow(BW1b);
pause(1);
title(FN1a, 'Interpreter', 'none');
bw3time = toc(bwstart)

% remove small areas again with a similar filtering size as
% before - but this should eliminate small background regions or other
% regions that are not cells

BW1b = bwareaopen(BW1b,InverseBWMaxPix);

figure('Numbertitle', 'off','Name','Function: InverseBW.m and bwareaopen');
imshow(BW1b);
title(FN1a, 'Interpreter', 'none');

clearvars -except BW1a BW1b
