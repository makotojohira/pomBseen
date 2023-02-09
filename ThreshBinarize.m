
% Apply Otsu's thresholding method, binarize, and return resulting BW
% image.

function [BW1b] = ThreshBinarize(R1a, FN1a, ThreshBinSensitivity, ThreshBinNeighbrhd, InverseBWMaxPix)

THD = adaptthresh(R1a, ThreshBinSensitivity, 'NeighborhoodSize', ThreshBinNeighbrhd);
% ThreshBinSensitivity = 0.55 defined in the Parameters function
% ThreshBinNeighbrhd = [15 15] defined in the Parameters function

BW1 = imbinarize(R1a, THD);
% So far this performs best on cells with side illumination improving with 
% neighborhoodsize of 21,21 -> 19,19 -> 17,17 -> 15,15 - improvememnt stops
% at 13,13 - so use 15,15



% Now I want to see what this looks like...

figure('Numbertitle', 'off','Name','Function: ThreshBinarize.m');
imshow(BW1);
pause(1);
title(FN1a, 'Interpreter', 'none');


% Since the halo is the key to segmenting I want to emphasize the halo's
% white line, per these suggestions: https://www.mathworks.com/matlabcentral/answers/479905-boost-enhance-white-pixels?s_tid=prof_contriblnk

BW1a = bwmorph(BW1, 'thicken', 1); % thicken the halo to help close small gaps
BW1a = bwareaopen(BW1a,InverseBWMaxPix); % eliminate small white specks within cells
BW1b = bwmorph(BW1a, 'bridge'); % bridge small gaps in the halo

% Let's visualize the effects of these morphological operations

figure('Numbertitle', 'off','Name','Function: ThreshBinarize.m + bwmorph');
imshow(BW1b);
pause(1);
title(FN1a, 'Interpreter', 'none');




clearvars -except BW1 R1a BW1b
