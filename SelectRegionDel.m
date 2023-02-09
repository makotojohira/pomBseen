
% Assuming the filtering was not perfect, manually select regions to delete
% and return image.

% Incorporate SegmentNum within this function to help with renumbering cells.
% Modify so user can enter multiple numbers separated by space.

function [BW1e,CC] = SelectRegionDel(BW1e,CC,FN1a)

% Ask the user for a region number plotted in the previous figure and 
% confirm that it is deleted.  Refer to the documentation on bwconncomp 
% here:
% https://www.mathworks.com/help/images/ref/bwconncomp.html

while(1);
titleBar1 = 'Select cell number';
userInput1 = inputdlg({'Please select cell number(s) to delete - put space between numbers'},titleBar1);
if isempty(userInput1),return,end;

N = cell2mat(userInput1)
N = str2num(N)
[m,n] = size(N)

for n = 1:n
    BW1e(CC.PixelIdxList{N(n)}) = 0;
end


%figure;
%imshow(BW1e);

CC      = bwconncomp(BW1e,4);
stats   = regionprops(CC,{'Area' 'Solidity'});
Area    = [stats.Area];
mask    = Area < 500 | Area > 100000;
CC.PixelIdxList(mask) = [];
CC.NumObjects   = length(CC.PixelIdxList);
Area(mask)      = [];
BW1e      = false(size(BW1e));
BW1e(vertcat(CC.PixelIdxList{:})) = true;

labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled,'spring','c','shuffle');

s = regionprops(CC,'centroid');
centroid = cat(1,s.Centroid);
%label = CC.PixelIdxList{1:end};
n = CC.NumObjects;
figure('Numbertitle', 'off','Name','Function: SegmentNum.m'); hold on
imshow(RGB_label);
pause(1);
%hold on;
%plot(centroid(:,1),centroid(:,2),'b*') % Confirm location of centroids
for n=1:n;
    text(centroid(n,1),centroid(n,2),sprintf('%d',n),'HorizontalAlignment','center');
end
title(FN1a, 'Interpreter', 'none');
hold off;



yn = menu('Do you want to enter another value for N?','Yes','No');
if yn==2 | yn==0;
    break
end

end

clearvars -except BW1e CC 

