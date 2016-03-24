clc;clear all;close all;
ckt = imread('scan1.jpg');
%ckt = imresize(ckt,[2500 2500]);
bwckt = im2bw(ckt);
% bwckt = ~imdilate(imrotate(~bwckt, 33, 'loose'),strel('disk',1));
% %simulate rotation
% imshow(bwckt)
% imshow(edge(im2bw(ckt), 'canny', .5, 6))
% imshow(bwckt);
%% Denoising
ratio = round(mean(sum(bwckt))/400)
denoise2 = wiener2(bwckt,[ratio ratio]);
imshow(imdilate(bwckt, strel('disk',round(ratio/6))));

%%
% test = imerode(bwckt,strel('disk',round(1.5*ratio)));
% imshow(test)
figure;
test2 = imdilate(bwckt, strel('disk',round(ratio/6)));
figure;
imshow(test2)
test3 = imerode(bwckt,strel('disk',ratio));
figure;
imshow(test3)
%%
sample_dilation_denoised = imdilate(imerode(bwckt,strel('disk',1)),strel('disk',1));
[sample_label, num_islands] = bwlabel(sample_dilation_denoised);

for aa = 1:num_islands
    if( sum(sum(sample_label == aa)) <= 125 )
        sample_label(sample_label == aa) = 0;
    end
end
bwckt_denoised = im2bw(sample_label);
figure; imshow(bwckt_denoised); title('Denoised');

bwckt_skel= bwmorph(imdilate(~bwckt_denoised,strel('disk',1)),'thin',Inf);
%bwckt_skel = bwmorph(~bwckt_denoised,'skel',Inf);
figure; imshow(bwckt_skel); title('Skeletonized');

BW = imdilate(bwckt_skel,strel('line',5,90));
%BW = edge(bwckt_skel, 'canny')

figure;imshow(BW);hold on; title('Dilated w/ Initial Lines');

%circle jerk detection
% [centers, radii, metric] = imfindcircles(BW,[20 40]);
% viscircles(centers, radii(1),'EdgeColor','b');

% Hough transform line detection
[H,theta,rho] = hough(BW, 'Theta', [-90, 82]);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW,theta,rho,P,'FillGap',3,'MinLength',floor(min(size(BW))/2.4));

lineVec = horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2));

%do the vertical lines (sooooo stupid)
BW2 = imdilate(bwckt_skel',strel('line',5,90));
[H,theta,rho] = hough(BW2, 'Theta', [-90, 82]);
P = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(BW2,theta,rho,P,'FillGap',3,'MinLength',40);

lineVec = vertcat(lineVec, horzcat(fliplr(vec2mat([lines.point1],2)),fliplr(vec2mat([lines.point2],2))));

for k = 1:size(lineVec, 1)
   xy = [lineVec(k,1:2); lineVec(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end
numLines = length(lineVec)

disp 'Joining nearby points'
allpts = vertcat(lineVec(:,1:2),lineVec(:,3:4));
ptsJoined = 0;
thresh = 15;
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

figure;imshow(~BW);hold on;
for k = 1:size(lines2, 1)
   xy = [lines2(k,1:2); lines2(k,3:4)];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

% remove all pts from image within 3px of lines
thresh = 5;
for xxx = 1:size(BW,2)
    for yyy = 1:size(BW,1)
        if(BW(yyy,xxx) == 1)
            for lx = 1:size(lines2,1)
               v = [lines2(lx,4)-lines2(lx,2), -(lines2(lx,3)-lines2(lx,1))];
               midpt = [lines2(lx,1)+lines2(lx,3), lines2(lx,2)+lines2(lx,4)]/2;
               dm = sqrt((midpt(1)-xxx)^2 + (midpt(2)-yyy)^2);
               len = sqrt((lines2(lx,1)-lines2(lx,3))^2 + (lines2(lx,2)-lines2(lx,4))^2);
               r = [lines2(lx,1)-xxx, lines2(lx,2)-yyy];
               dist = abs(dot(v/norm(v),r));
               if(dist<thresh && dm <= len/2)
                   BW(yyy,xxx) = 0;
               end
            end
        end
    end
end
% imshow(BW)
%%
[foo,bar]=bwlabel(BW);

cent = regionprops(foo,'Centroid');
cent2 = reshape([cent.Centroid],2,[])';
thresh = 50;
countThresh = 20;
for ii = 1:length(cent2)
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

islands = regionprops(foo, 'BoundingBox');
islandCount = size(islands, 1)

components = zeros(1,islandCount);

padding = 10;
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
    wx = min(size(bwckt,2)-x1, rects(3)+2*padding);
    wy = min(size(bwckt,1)-y1, rects(4)+2*padding);
    
    components = [components imcrop(bwckt, [xx1 yy1 wx wy])];
%     components = [components imcrop(bwckt, islands(k).BoundingBox)];
end

figure;

for ii = 1:length(components)
   subplot(3, ceil(islandCount/3), ii);
   imshow(components{ii});
end

%%

%corner detection
% c = corner(BW, 'Harris');
% plot(c(:,1),c(:,2),'b*')
% keptCorners = [];
% for m = 1:length(c)
%    cx = c(m,1);
%    cy = c(m,2);
%    corn = [c(m,1), c(m,2)];
%    for k = 1:length(lines)
%         distance = abs(det([corn-lines(k).point1;lines(k).point1-lines(k).point2]))/norm(lines(k).point1-lines(k).point2);
%         if distance < 1
%             distance
%             keptCorners = [keptCorners;c(m,:)];
%         end
%    end 
% end
% 
% hold on;
% plot(keptCorners(:,1),keptCorners(:,2),'m.')


