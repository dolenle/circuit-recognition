function [ components_preprocess, pins , netlist_mtx] = segment3( ckt_file, min_line)
%%
%ckt_file = 'Scan 6.jpg'; min_line = 200;
clc;
% Read In Scanned Image & Crop Border
scaned_image = imread(ckt_file);
cropped_scan = imcrop(scaned_image, [10 10 size(scaned_image,2)-20 size(scaned_image,1)-20]); %crop 10px edges
    %figure;imshow(cropped_scan); title('cropped_scan');
% Resize scan image and invert image
resized_scan = imresize(cropped_scan,[4096 4096]);
resized_scan = ~im2bw(resized_scan,0.6);
    figure; imshow(resized_scan); title('Resized Inverted Scanned Image')
    hold on;
% Thin Scanned Image 
dilate_resized_scan = imdilate(resized_scan, strel('disk', 10));
thin_scan = bwmorph(dilate_resized_scan, 'thin', Inf);
    %figure;imshow(thin_scan); title('thin_scan');
% Denoising by removing island size
    % Remove islands smaller than 1/10th of average
    %avgIslandSize = mean(cellfun(@length, CC.PixelIdxList))
removed_specks = 0;
CC = bwconncomp(thin_scan);
denoised_scan = thin_scan;
for aa = 1:CC.NumObjects
    %if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
    if(length(CC.PixelIdxList{aa}) < 500)
        denoised_scan(CC.PixelIdxList{aa}) = 0;
        removed_specks = removed_specks+1;
    end
end
less_denoised_scan = thin_scan;
for aa = 1:CC.NumObjects
    %if(length(CC.PixelIdxList{aa}) < avgIslandSize/10)
    if(length(CC.PixelIdxList{aa}) < 50)
        less_denoised_scan(CC.PixelIdxList{aa}) = 0;
    end
end
removed_specks;
    %figure;imshow(denoised_scan); title('denoised scan');
    %figure;imshow(less_denoised_scan); title('less denoised scan');

% Thicken Images and Transpose for Hough Transform
thick_denoised_scan = imdilate(denoised_scan, strel('rectangle', [1,20])); % Used 'disk' @ 930PM
thick_denoised_scan_rot = imdilate(denoised_scan', strel('rectangle', [1,20]));
%thick_denoised_scan = thick_denoised_scan_rot;
    %figure; imshow(thick_denoised_scan); title('thick_denoised_scan');

% Compute Hough Transform for Line Detection
[H,theta,rho] = hough(thick_denoised_scan);
[H_rot,theta_rot,rho_rot] = hough(thick_denoised_scan_rot);
% Visualize Hough Transform
    % figure;imshow(imadjust(mat2gray(H)),'XData',theta,'YDATA',rho,'InitialMagnification','fit');
    % axis on;axis normal,hold on; colormap(hot);
peaks = houghpeaks(H,100,'threshold',ceil(0.3*max(H(:))));
peaks_rot = houghpeaks(H_rot,100,'threshold',ceil(0.3*max(H_rot(:))));
    % plot(theta(peaks(:,2)),rho(peaks(:,1)),'s','color','black');

% Find Lines of Hough Transform
lines = houghlines(thick_denoised_scan,theta,rho,peaks,'FillGap',3,'MinLength',min_line);
pre_deleted_lines = lines;
lines_rot = houghlines(thick_denoised_scan_rot,theta_rot,rho_rot,peaks_rot,'FillGap',3,'MinLength',min_line);
pre_deleted_lines_rot = lines_rot;
    %hold off;

%figure; imshow(~imdilate(denoised_scan, strel('disk', 5))); title('not thick denoised scan'); hold on;
for ii = 1:2
    if (ii == 2)
        lines = lines_rot;
    end
    % Remove angled lines (i.e. BJT) and any lines which intersect them=
    angles = [lines.theta];
    slantInd = find(angles < 50 & angles > 40 | angles > -50 & angles < -40);
    slantLines = horzcat(vec2mat([lines(slantInd).point1],2),vec2mat([lines(slantInd).point2],2)); %End points of angled lines
    row_index = [];
    output_intersections = lineSegmentIntersect(horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2)), slantLines);
    [row_index,col_index] = find(output_intersections.intAdjacencyMatrix);
    
    slantInd = unique([slantInd row_index']);
    % Delete Slanted & Intersected Lines
    lines(slantInd) = [];
    linesRemoved = length(slantInd);
    lineVec = horzcat(vec2mat([lines.point1],2),vec2mat([lines.point2],2));
    %disp 'Joining nearby points'
    allpts = vertcat(lineVec(:,1:2),lineVec(:,3:4)); % everything x two
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
    %ptsJoined
    
    % Remove duplicate lines
    %disp 'Removing duplicate lines'
    line_mtx_no_duplicates = unique(horzcat(allpts(1:size(allpts,1)/2,:),allpts(size(allpts,1)/2+1:end,:)),'rows');
    %linesRemoved = size(allpts,1)/2-size(line_mtx_no_duplicates,1)
    if (ii == 1)
        line_mtx = line_mtx_no_duplicates;
    else
        swaped_xy = zeros(size(line_mtx_no_duplicates));
        swaped_xy(:,1) = line_mtx_no_duplicates(:,2);
        swaped_xy(:,2) = line_mtx_no_duplicates(:,1);
        swaped_xy(:,3) = line_mtx_no_duplicates(:,4);
        swaped_xy(:,4) = line_mtx_no_duplicates(:,3);
        line_mtx = [line_mtx; swaped_xy];
    end
end

for k = 1:size(line_mtx, 1)
    xy = [line_mtx(k,1:2); line_mtx(k,3:4)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    % Plot beginnings and ends of lines
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

% Erasing Lines From Original Image
%ckt_no_line = denoised_scan;
ckt_no_line = less_denoised_scan;
%ckt_no_line = imdilate(denoised_scan, strel('disk', 5));
point_mtx = line_mtx;
erase_thickness = 30;
erase_limit = 10;
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

    for xx = (min(x1,x2)+erase_limit):((max(x1,x2))-erase_limit)
        yy = round(mm*xx+bb);
        yy_min = yy-erase_thickness;
        yy_max = yy+erase_thickness;
        if (horiz == 1)
            ckt_no_line([yy_min:yy_max],xx)=0;
        else
            ckt_no_line(xx,[yy_min:yy_max])=0;
        end
    end
end
%figure; imshow(ckt_no_line);
%figure; imshow(imdilate(ckt_no_line, strel('disk', 10)));

% Making The Erased Lines Touch
ckt_no_line = imdilate(ckt_no_line, strel('disk', 10));
ckt_no_line = bwmorph(ckt_no_line, 'thin', Inf);

% Extracting Components and Nodes
[label,number_of_islands]=bwlabel(ckt_no_line);
nodes = label;
for aa = 1:number_of_islands
            if( sum(sum(label == aa)) <= 49 )
                label(label == aa) = 0;
            else
                nodes(nodes == aa) = 0;
            end
end
%figure;imshow(imdilate(nodes, strel('disk', 5)));
%figure;imshow(imdilate(label, strel('disk', 5)));

% Combine Nearby Sections
cent = regionprops(label,'Centroid');
cent = reshape([cent.Centroid],2,[])';
%hold on;scatter(cent(:,1),cent(:,2));hold off;

thresh = 400;
countThresh = 100;
combined_labels = label;
for ii = 1:length(cent)
    x = cent(ii,:);
    y = cent(ii+1:end,:);
    distances = sqrt((y(:,1)-x(1)).^2 + (y(:,2)-x(2)).^2);
    nearby = ii+find(distances < thresh);
    for jj = 1:length(nearby)
       combined_labels(combined_labels == nearby(jj)) = ii; 
    end
end
ii = 1;
while 1==1
    if ii > max(max(combined_labels))
        break
    end
    if sum(combined_labels(:) == ii) < countThresh 
        combined_labels(combined_labels == ii) = 0;
        combined_labels(combined_labels > ii) = combined_labels(combined_labels > ii)-1;
        ii = ii-1;
    end
    ii = ii+1;
end
hold on;
islands = regionprops(combined_labels, 'BoundingBox');
num_of_components = size(islands, 1)

components = zeros(1,num_of_components);

padding = 100;
rectCollection = [];
components = {};
only_wires = denoised_scan;
top = zeros(num_of_components,2); bot = zeros(num_of_components,2); 
left = zeros(num_of_components,2); right = zeros(num_of_components,2);
for k = 1:num_of_components
    rects = islands(k).BoundingBox;
    x1 = rects(1);
    y1 = rects(2);
    x2 = x1 + rects(3);
    y2 = y1 + rects(4);
    x = [x1 x2 x2 x1 x1];
    y = [y1 y1 y2 y2 y1];
    plot(x, y, 'LineWidth', 2);
    xx1 = max(0, rects(1)-padding);
    yy1 = max(0, rects(2)-padding);
    wx = min(size(resized_scan,2)-x1, rects(3)+2*padding);
    wy = min(size(resized_scan,1)-y1, rects(4)+2*padding);

    % Roundapalooza
    xx1=round(xx1);yy1=round(yy1);wx=round(wx);wy=round(wy);
    
    components = [components imcrop(resized_scan, [xx1 yy1 wx wy])];
    only_wires(yy1:yy1+wy,xx1:xx1+wx) = 0;
    
    if (~isempty(find(denoised_scan(yy1,xx1:xx1+wx))))
        top(k,:) = [xx1+find(denoised_scan(yy1,xx1:xx1+wx)),yy1];
    end
    if (~isempty(find(denoised_scan(yy1+wy,xx1:xx1+wx))))
        bot(k,:) = [xx1+find(denoised_scan(yy1+wy,xx1:xx1+wx)),yy1+wy];
    end
    if (~isempty(find(denoised_scan(yy1:yy1+wy,xx1))))
        left(k,:) = [xx1,yy1+find(denoised_scan(yy1:yy1+wy,xx1))];
    end
    if (~isempty(find(denoised_scan(yy1:yy1+wy,xx1+wx))))
        right(k,:) = [xx1+wx,yy1+find(denoised_scan(yy1:yy1+wy,xx1+wx))];
    end

%     pins(ii) = top+bot+left+right;

    
end

%figure;imshow(denoised_scan);title('Original Thinned Ckt');
%figure;imshow(imdilate(only_wires,strel('disk',5)));
% Get Wire Islands
labeled_wires = bwlabel(only_wires);
end_points = bwmorph(only_wires,'endpoints');
labeled_end_points = double(end_points);
labeled_end_points(end_points==1) = labeled_wires(end_points==1);
[epY,epX, ep_values] = find(labeled_end_points);

%figure; imshow(~only_wires); hold on;
for ii = 1:length(ep_values)
    scatter(epX(ii),epY(ii),[],'filled');
    tt = text(epX(ii)+20,epY(ii)-20,cellstr(num2str(ep_values(ii))));
    set(tt,'Color',[1 0 0]);
    set(tt,'FontSize',16);
end

% Individual Componnets Denosied
% hold off;
figure; title('Extracted Components');
components_preprocess = cell(1);
for ii = 1:length(components)
   CC_components = bwconncomp(components{ii});
   for aa = 1:CC_components.NumObjects
        if(length(CC_components.PixelIdxList{aa}) < 50)
            components{ii}(CC_components.PixelIdxList{aa}) = 0;
        end
   end
    subplot(3, ceil(num_of_components/3), ii);
    imshow(imresize(components{ii},[250 250]));
    components_preprocess{1} = cat(3,components_preprocess{1}, imresize(components{ii},[250 250]));
end
pins = zeros(1,size(components_preprocess{1},3));
%%
% Extracting Number of Pins Per Component
for ii = 1:size(components_preprocess{1},3)
    top_flag = (0<sum(components_preprocess{1}(1,:,ii)));
    bot_flag = (0<sum(components_preprocess{1}(end,:,ii)));
    left_flag = (0<sum(components_preprocess{1}(:,1,ii)));
    right_flag = (0<sum(components_preprocess{1}(:,end,ii)));
    pins(ii) = top_flag+bot_flag+left_flag+right_flag;
end

% Componennt Centriods? Lolol
component_centriods = regionprops(combined_labels,'Centroid');
component_centriods = reshape([component_centriods.Centroid],2,[])';
%figure; hold on; imshow(~resized_scan); hold on;
% for ii = 1:length(component_centriods)
%     scatter(component_centriods(ii,1),component_centriods(ii,2),[],'filled');
%     text(component_centriods(ii,1),component_centriods(ii,2),cellstr(num2str(ii)));
% end


% Generating netlist matrix
netlist_mtx = zeros(num_of_components, 5);
end_points_xy = [epX epY];
nets = ep_values;
% connect ground to zero
for ii = 1:num_of_components
    if(pins(ii)==1)
        point = [top(ii,:);bot(ii,:);left(ii,:);right(ii,:)];
        point = point(any(point,2),:);
        oldLabel = nets(find(abs(sum(bsxfun(@minus, end_points_xy, point),2))<5));
        nets(nets == oldLabel) = 0;
    end
end
% assign nodes
for ii = 1:num_of_components
    if (~(pins(ii) == 1))
        if(pins(ii) == 2) %2 pin devices
            if(right(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, right(ii,:)),2));
                netlist_mtx(ii,3) = nets(find(dist<5));
            end
            if(left(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, left(ii,:)),2));
                netlist_mtx(ii,4) = nets(find(dist<3));
            end
            if(top(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, top(ii,:)),2));
                netlist_mtx(ii,3) = nets(find(dist<3));
            end
            if(bot(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, bot(ii,:)),2));
                netlist_mtx(ii,4) = nets(find(dist<3));
            end
        elseif(pins(ii) == 3) %3-pin devices (xstr)
            if(right(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, right(ii,:)),2));
                netlist_mtx(ii,4) = nets(find(dist<3));
            end
            if(left(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, left(ii,:)),2));
                netlist_mtx(ii,4) = nets(find(dist<3));
            end
            if(top(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, top(ii,:)),2));
                netlist_mtx(ii,3) = nets(find(dist<3));
            end
            if(bot(ii,1))
                dist = abs(sum(bsxfun(@minus, end_points_xy, bot(ii,:)),2));
                netlist_mtx(ii,5) = nets(find(dist<3));
            end
        end
    end
end
%%
end