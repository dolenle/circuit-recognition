% clc;close all;

function [Xind, Yind] = getGrid(expected_rows, expected_cols, cktImg)
% disp(path);
% figure; imshow(cktImg);
% scan = ~im2bw(imread(path),0.95);
enableDebug = true;
cktImg = ~cktImg;
colPeaks = sum(cktImg);
[Xpeaks, Xind] = findpeaks(colPeaks, 'minpeakheight', max(colPeaks)*0.7);
rowPeaks = sum(cktImg, 2);
sortPeaks = unique(rowPeaks);
[Ypeaks, Yind] = findpeaks(rowPeaks, 'minpeakheight', sortPeaks(end-1)/5,'sortstr','descend');

disp 'verifying grid'
verifyMargin = 12;
spaceMargin = 16;
Xrange = [Xind-4 Xind-3 Xind-2 Xind-1 Xind Xind+1 Xind+2 Xind+3 Xind+4];

badY = [];
for ii = 1:length(Yind)
    row = Yind(ii);
    above = row-verifyMargin;
    below = row+verifyMargin;
    aboveS = row-spaceMargin;
    belowS = row+spaceMargin;
    if(above > 0 && below < size(cktImg,1))
        if((sum(cktImg(above,Xrange)) + sum(cktImg(below,Xrange)) < 4) ...
            || (sum(cktImg(aboveS,setdiff(1:size(cktImg,2),Xrange)))>80 || sum(cktImg(belowS,setdiff(1:size(cktImg,2),Xrange)))>80))
            disp 'invalid row at'; disp(Yind(ii))
            badY = [badY ii];
        end
    end
end
Yind(badY) = [];
Ypeaks(badY) = [];

Yrange = [Yind-4 Yind-3 Yind-2 Yind-1 Yind Yind+1 Yind+2 Yind+3 Yind+4];
badX = [];
for ii = 1:length(Xind)
    col = Xind(ii);
    left = col-verifyMargin;
    right = col+verifyMargin;
    leftS = col-spaceMargin;
    rightS = col+spaceMargin;
    if(left > 0 && right < size(cktImg,2))
        if((sum(cktImg(Yrange,left)) + sum(cktImg(Yrange,right)) < 4) ...
            || (sum(cktImg(setdiff(1:size(cktImg,1),Yrange),leftS))>80 || sum(cktImg(setdiff(1:size(cktImg,1),Yrange),rightS))>80))
            disp 'invalid column at'; disp(Xind(ii))
            badX = [badX ii];
        end
    end
end
Xind(badX) = [];
Xpeaks(badX) = [];
%%
disp 'combining nearby points'
Xrange = 15;
Xremove = [];
for ii = 1:length(Xind)
    nearby = find(Xind(ii)-Xrange < Xind & Xind(ii)+Xrange > Xind);
    if(length(nearby) > 1)
        Xremove = [Xremove nearby(2:end)];
        Xpeaks(nearby) = max(Xpeaks(nearby));
    end
end
Xind(unique(Xremove)) = [];
Xpeaks(unique(Xremove)) = [];

Yrange = 15;
Yremove = [];
for ii = 1:length(Yind)
    nearby = find(Yind(ii)-Yrange < Yind & Yind(ii)+Yrange > Yind)';
    if(length(nearby) > 1)
        Yremove = [Yremove nearby(2:end)];
        Ypeaks(nearby) = max(Ypeaks(nearby));
    end
end
Yind(unique(Yremove)) = [];
Ypeaks(unique(Yremove)) = [];

if(length(Yind) ~= expected_rows)
    figure; imshow(cktImg);
    figure; hold on;
    plot(linspace(1,size(cktImg,2),size(cktImg,2)), colPeaks);
    plot(Xind, Xpeaks, 'r*');
    figure; hold on;
    plot(linspace(1,size(cktImg,1),size(cktImg,1)), rowPeaks);
    plot(Yind, Ypeaks, 'go');
    view(-90, 90);
    set(gca,'Xdir','reverse','Ydir','reverse')
    fprintf(2,'ERROR: Incorrect number of rows\n')
    if(enableDebug)
        disp 'STOP';
    end
end
if(length(Xind) ~= expected_cols)
    figure; imshow(cktImg);
    figure; hold on;
    plot(linspace(1,size(cktImg,2),size(cktImg,2)), colPeaks);
    plot(Xind, Xpeaks, 'r*');
    figure; hold on;
    plot(linspace(1,size(cktImg,1),size(cktImg,1)), rowPeaks);
    plot(Yind, Ypeaks, 'go');
    view(-90, 90);
    set(gca,'Xdir','reverse','Ydir','reverse')
    fprintf(2,'ERROR: Incorrect number of columns\n')
    if(enableDebug)
        disp 'STOP';
    end
end
%% Plot
% figure; hold on;
% plot(linspace(1,size(cktImg,2),size(cktImg,2)), colPeaks);
% plot(Xind, Xpeaks, 'r*');
% figure; hold on;
% plot(linspace(1,size(cktImg,1),size(cktImg,1)), rowPeaks);
% plot(Yind, Ypeaks, 'go');
% view(-90, 90);
% set(gca,'Xdir','reverse','Ydir','reverse')
%%
% width = 10;
% figure; imshow(cktImg);
% columns = Xind;
% rows = Yind;
% % Deleting Rows and Columns
% for kk = -width:width
%     cktImg(rows+kk,:) = 1;
%     cktImg(:,columns+kk) = 1;
% end
% figure; imshow(cktImg);
%% Testing
% row = 3250;
% verifyMargin = 12;
% spaceMargin = 16;
% above = row-verifyMargin;
% below = row+verifyMargin;
% aboveS = row-spaceMargin;
% belowS = row+spaceMargin;
% sum(cktImg(above,Xrange)) + sum(cktImg(below,Xrange))
% sum(cktImg(aboveS,setdiff(1:size(cktImg,2),Xrange)))
% sum(cktImg(belowS,setdiff(1:size(cktImg,2),Xrange)))
