clc;clear all;close all;
ckt = imread('scan1a.jpg');
%ckt = imresize(ckt,[2500 2500]);
bwckt = im2bw(ckt);
% bwckt = ~imdilate(imrotate(~bwckt, 33, 'loose'),strel('disk',1));
% %simulate rotation
% imshow(bwckt)
% imshow(edge(im2bw(ckt), 'canny', .5, 6))
imshow(bwckt);
imshow(imerode(bwckt, strel('disk',2)));
%% Denoising
ratio = round(mean(sum(bwckt))/200); %estimate size

%remove islands smaller than 1/10th of average
removed_specks = 0;
CC = bwconncomp(~bwckt);
avgIslandSize = mean(cellfun(@length, CC.PixelIdxList))
for aa = 1:CC.NumObjects
    if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
        bwckt(CC.PixelIdxList{aa}) = 1;
        removed_specks = removed_specks+1;
    end
end
removed_specks
imshow(bwckt)

thickened = imerode(bwckt, strel('disk',round(ratio)));
figure; imshow(thickened);

%%
edges = edge(thickened);
corn = corner(thickened, 5);
[H,theta,rho] = hough(edges);
peaks = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
lines = houghlines(edges,theta,rho,peaks);
figure, imshow(bwckt), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end

plot(corn(:,1), corn(:,2), 'r*');