
% Delete cells contacting image border, repeat filtering, and return image.

function [BW1c] = ClearBorder(BW1b,FN1a,ClearBorderConn,ClearBorderMaxPix)

% Remove objects touching the border

BW1c = imclearborder(BW1b,ClearBorderConn);

% Variable ClearBorderConn = 8 defined in the Parameters function

% Repeat the filtering of small bright spots within cells - but now 
% excessive noise outside of cells.  I can filter sizes larger now so
% increase from 600 pixels in previous block... try double

BW1c = bwareaopen(BW1c,ClearBorderMaxPix);
% ClearBorderMaxPix = 1200 defined in the Parameters function

% Display results of ClearBorder function:

figure('Numbertitle', 'off','Name','Function: ClearBorder.m');
imshow(BW1c);
pause(1);
title(FN1a, 'Interpreter', 'none');

clearvars -except BW1c
