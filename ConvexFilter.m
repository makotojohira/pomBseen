
% Use convex polygon to filter out non-cell regions that display
% significant concavity, and return image.  Also use aspect ratio as an
% important modifier to convexity in order to retain long cells which
% happen to be bent due to their length.

% Modify the filtering using a linear combination of the area to convex 
% area ratio and the aspect ratio

function [CCstats,BW1e] = ConvexFilter(CC,BW1d,FN1a,ConvexFilterSlope,ConvexFilterIntercept)

CCstats = regionprops('table', CC, 'Area', 'ConvexArea', 'MajorAxisLength', 'MinorAxisLength');
[r,c] = size(CCstats);
CCstats.idx = (1:r)';
CCstats.AreaRatio = CCstats.ConvexArea ./ CCstats.Area;
CCstats.AspectRatio = CCstats.MajorAxisLength ./ CCstats.MinorAxisLength;

figure;
plot(CCstats.AreaRatio, CCstats.AspectRatio, 'g.');
xlabel('Area Ratio (Convex Area / Area)');
ylabel('Aspect Ratio (Major Axis L / Minor Axis L)');
hold on;
s = regionprops(CC,'centroid');
centroid = cat(1,s.Centroid);

for n=1:r
    text(CCstats.AreaRatio(n), CCstats.AspectRatio(n),sprintf('%d',n),'HorizontalAlignment','right');
end
drawnow;
title(FN1a, 'Interpreter', 'none');

% Now plot line of equation from below

MinAR = min(CCstats.AreaRatio);
MaxAR = max(CCstats.AreaRatio);
x = MinAR:MaxAR;
y = (ConvexFilterSlope * x) - ConvexFilterIntercept;
plot(x,y);

hold off;

% Based on the data shown in these plots, the line that passes through the
% points (1.05,1) and (1.4,5.5) appears to maximize the number of cells
% and minimize the artifacts retained above the line, while maximizing the
% artifacts and minimizing the cells discarded below the line.  The slope
% of the line is 4.5/0.35 = 12.8571.  The equation is y = 12.8571x + b; to
% solve for b substitute one of the points: 1 = (12.8571)(1.05) + b so we
% have b = -12.5.  Therefore:   y = 12.8571x - 12.5

% In our example, x = area ratio and y = aspect ratio.
% We can make this into a logical function such that if a region's aspect
% ratio > y for a given area ratio, that indexed region has a value of 1,
% otherwise it has a value of 0.

% Like the previous version of ConvexFilter, I'll implement this in a table
% or matrix instead of as a logical filter.

CCstats.y = (ConvexFilterSlope * CCstats.AreaRatio) - ConvexFilterIntercept;
%ConvexFilterSlope = 12.8571 defined in the Parameters function
%ConvexFilterIntercept = 12.5 defined in the Parameters function

% Now I want to make the pixels of an index region have a value = 1 IF the
% aspect ratio > y ELSE value = 0

mask = CCstats.AspectRatio < CCstats.y;
CC.PixelIdxList(mask) = [];
CC.NumObjects   = length(CC.PixelIdxList);

labeled = labelmatrix(CC);
RGB_label = label2rgb(labeled,'spring','c','shuffle');
figure('Numbertitle', 'off','Name','Function: Segment.m');
imshow(RGB_label);
pause(1);
title(FN1a, 'Interpreter', 'none');


BW1e = false(size(BW1d));
BW1e(vertcat(CC.PixelIdxList{:})) = true;

figure('Numbertitle', 'off','Name','Function: ConvexFilter.m');
imshow(BW1e);
title(FN1a, 'Interpreter', 'none');
drawnow;

clearvars -except CCstats BW1e

%%
%CCstats.PixelValue = double(CCstats.AspectRatio > CCstats.y)
%zeroIdx = CCstats.idx(CCstats.PixelValue == 0)
%numZero = size(zeroIdx)
%for n = 1:numZero
%    BW1e = BW1e
%    N = zeroIdx(n)
%    BW1e(CC.PixelIdxList{N}) = 0
%end



%%
%CCstats.idx = (1:end)
%CCstats.AreaRatio = CCstats.Area(1:end) / CCstats.ConvexArea(1:end)
%CCstats.AspectRatio = CCstats.MajorAxisLength(1:end) / CCstats.MinorAxisLength(1:end)

%area = regionprops(CC, 'Area');
%convexArea = regionprops(CC, 'ConvexArea');
%[r,c] = size(area);
%areaM = cell2mat(struct2cell(area)); % Convert structure array to matrix;
%convAM = cell2mat(struct2cell(convexArea));
%areaM = areaM'; % Convert rows to colums
%convAM = convAM';
%areasM = [areaM convAM]; % Combine into single matrix
%areasM(:,3) = [1:r]; % Add a 3rd column with index numbers
% Create 4th column with a threshold of 85% of convex area - determined by
% trial and error with one image so far... 75% seems to work OK with
% longest cells but also removes significantly bent cells and fails to
% remove regions that are mostly straight-sided polygons
%areasM(:,4) = areasM(:,2) .* 0.85;
% Now I want to find the index number of cells in which the area is less
% than the threshold
%whichRows = areasM(areasM(:,1) < areasM(:,4),:);
%delIdx = whichRows(:,3);
% Now I want to apply the index to the image and delete those cells whose
% area is less than the 85% convex area threshold
%BW1e = BW1d;
%rIdx = size(delIdx);
%for n = 1:rIdx;
%    BW1e = BW1e;
%    N = delIdx(n);
%    BW1e(CC.PixelIdxList{N}) =0;
%end
%figure('Numbertitle', 'off','Name','Function: ConvexFilter.m');
%imshow(BW1e);
%title(tiffFilename, 'Interpreter', 'none');
%drawnow;
% 85% threshold works great with smallest cells, not so good with larger
% cells, especially if they are bent as longer cells typically are.

