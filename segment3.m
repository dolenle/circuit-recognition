clc;clear all;close all;
ckt = imread('scans/scan3.jpg');
ckt = imcrop(ckt, [10 10 size(ckt,2)-20 size(ckt,1)-20]); %crop 10px edges
imshow(ckt);
ckt = imresize(ckt,[4096 4096]);
ckt = ~im2bw(ckt,0.6);
cktOrig = ckt;

ckt = imdilate(ckt, strel('disk', 10));
ckt = bwmorph(ckt, 'thin', Inf);

%remove islands smaller than 1/10th of average
removed_specks = 0;
CC = bwconncomp(ckt);
avgIslandSize = mean(cellfun(@length, CC.PixelIdxList))
for aa = 1:CC.NumObjects
    if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
        ckt(CC.PixelIdxList{aa}) = 0;
        removed_specks = removed_specks+1;
    end
end
removed_specks

ckt = imdilate(ckt, strel('disk', 10));
ckt_rot = ckt';

[H,theta,rho] = hough(ckt);
peaks = houghpeaks(H,100,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(ckt,theta,rho,peaks,'FillGap',3,'MinLength',avgIslandSize/10);
imshow(~ckt), hold on

%remove angled lines (i.e. BJT) and any lines which intersect them
disp 'remove slanted lines'
angles = [lines.theta];
slantInd = find(angles < 60 & angles > 30 | angles > -60 & angles < -30);
slantLines = horzcat(vec2mat([lines(slantInd).point1],2),vec2mat([lines(slantInd).point2],2));
INT = lineSegmentIntersect(horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2)), slantLines);
[r,c] = find(INT.intAdjacencyMatrix);
slantInd = unique([slantInd r']);
lines(slantInd) = [];
linesRemoved = length(slantInd);

lineVec = horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2));
disp 'Joining nearby points'
allpts = vertcat(lineVec(:,1:2),lineVec(:,3:4));
ptsJoined = 0;
thresh = 25;
for k = 1:size(allpts, 1)
    x = allpts(k,:);
    y = vertcat(allpts(1:k-1,:),allpts(k+1:end,:));
    distances = sqrt((y(:,1)-x(1)).^2 + (y(:,2)-x(2)).^2);
    pts = vertcat(x, y(find(distances < thresh),:));
    if(size(pts,1) > 1)
        pt = round(mean(pts));
        idx = find(distances < thresh);
        idx = horzcat(idx(idx<k)', k, idx(idx>k)');
        idx(idx>k) = idx(idx>k) + 1;
        allpts(idx,:) = repmat(pt,length(idx),1);
        ptsJoined = ptsJoined + length(pts)-1;
    end
end
ptsJoined

%remove duplicate lines
disp 'Removing duplicate lines'
temp = size(allpts,1)/2;
lines2 = unique(horzcat(allpts(1:size(allpts,1)/2,:),allpts(size(allpts,1)/2+1:end,:)),'rows');
linesRemoved = temp-size(lines2,1)

for k = 1:size(lines2, 1)
   xy = [lines2(k,1:2); lines2(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

%%
image = ckt;
point_mtx = lines2;
thickness = 30;
for ii = 1:size(point_mtx,1)
    if (abs(point_mtx(ii,1)-point_mtx(ii,3))< abs(point_mtx(ii,2)-point_mtx(ii,4)))
        x1 = point_mtx(ii,2); y1 = point_mtx(ii,1);
        x2 = point_mtx(ii,4); y2 = point_mtx(ii,3);
        horiz = 0;
    else
        x1 = point_mtx(ii,1);y1 = point_mtx(ii,2);
        x2 = point_mtx(ii,3);y2 = point_mtx(ii,4);
        horiz = 1;
    end
    mm = (y2-y1)/(x2-x1);
    bb = y1 - mm*x1;
    for xx = min(x1,x2):max(x1,x2)
        yy = round(mm*xx+bb);
        yy_min = yy-thickness;
        yy_max = yy+thickness;
        if (horiz == 1)
            image([yy_min:yy_max],xx)=0;
        else
            image(xx,[yy_min:yy_max])=0;
        end
    end
end
ckt = image;
%%
figure; imshow(~cktOrig);
[foo,bar]=bwlabel(ckt);

cent = regionprops(foo,'Centroid');
cent2 = reshape([cent.Centroid],2,[])';
thresh = 400;
countThresh = 100;
for ii = 1:length(cent2) %combine nearby sections
    x = cent2(ii,:);
    y = cent2(ii+1:end,:);
    distances = sqrt((y(:,1)-x(1)).^2 + (y(:,2)-x(2)).^2);
    nearby = ii+find(distances < thresh);
    for jj = 1:length(nearby)
       foo(foo == nearby(jj)) = ii; 
    end
end
ii = 1;
while 1==1
    if ii > max(max(foo))
        break
    end
    if sum(foo(:) == ii) < countThresh 
        foo(foo == ii) = 0;
        foo(foo > ii) = foo(foo > ii)-1;
        ii = ii-1;
    end
    ii = ii+1;
end
% for ii = 1:length(cent2)
%     if sum(foo(:) == ii) < countThresh 
%         kk = ii
%         foo(foo == ii) = 0;
%         foo(foo > ii) = foo(foo > ii)-1;
%         ii = ii-1
%     end
% end
hold on;
islands = regionprops(foo, 'BoundingBox');
islandCount = size(islands, 1)

components = zeros(1,islandCount);

padding = 30;
rectCollection = [];
components = {};
for k = 1:islandCount
    rects = islands(k).BoundingBox;
    x1 = rects(1);
    y1 = rects(2);
    x2 = x1 + rects(3);
    y2 = y1 + rects(4);
    x = [x1 x2 x2 x1 x1];
    y = [y1 y1 y2 y2 y1];
    plot(x, y, 'LineWidth', 2);
%     rectCollection(k,:,:,:) = [x1; y1; x2; y2];
%     test = foo;
%     test(foo~=k) = 0;
%     test2 = regionprops(test, 'Image');
%     components = [components test2.Image];
    xx1 = max(0, rects(1)-padding);
    yy1 = max(0, rects(2)-padding);
    wx = min(size(cktOrig,2)-x1, rects(3)+2*padding);
    wy = min(size(cktOrig,1)-y1, rects(4)+2*padding);
    
    components = [components imcrop(cktOrig, [xx1 yy1 wx wy])];
%     components = [components imcrop(bwckt, islands(k).BoundingBox)];
end

figure;

for ii = 1:length(components)
   subplot(3, ceil(islandCount/3), ii);
   imshow(components{ii});
end

%% Branch points

% skel = bwmorph(ckt, 'skel', Inf);
% figure; imshow(skel); hold on;
% ends = bwmorph(skel, 'endpoints');
% [row, column] = find(ends);
% endPts = [row column];
% branches = bwmorph(skel, 'branchpoints');
% [row, column] = find(branches);
% branchPts = [row column];
% plot(endPts(:,2), endPts(:,1), 'r*');
% plot(branchPts(:,2), branchPts(:,1), 'go');
% corn = corner(ckt, 50);
% plot(corn(:,1), corn(:,2), 'r*');