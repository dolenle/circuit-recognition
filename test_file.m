% Testing-File
% training_data_dir = dir('C:\Users\Clarakins\Documents\MATLAB\Senior Project\training_data');
% subfolders = training_data_dir(3:end);
% raw_images = cell(1,length(subfolders));
% for ii = 1:length(subfolders)
%     subfolder_dir = dir(['C:\Users\Clarakins\Documents\MATLAB\Senior Project\training_data\' subfolders(ii).name]);
%     subfolder_dir = subfolder_dir(3:end);
%     raw_images{ii} = cell(0,0);
%     for jj = 1:length(subfolder_dir)
%         raw_images{ii} = [raw_images{ii} imread(subfolder_dir(jj).name)];
%     end
% end
clc;clear all;close all;
load('training_images.mat')

test = training_images{1}{1};
figure; imshow(test)

test_denoised = imdilate(imerode(test,strel('disk',2)),strel('disk',2));
figure; imshow(test_denoised)

[test_labeled num]= bwlabel(test_denoised);
figure; imshow(test_labeled)
num

properties = {'Area','MajorAxisLength','MinorAxisLength','Centroid', ...
                'Orientation','Eccentricity','ConvexArea','FilledArea', ...
                'EquivDiameter','Extent','Solidity'};

features = [];
for ii = 1:length(properties)
    rp = regionprops(test_labeled,properties{ii});
    features = [features struct2array(rp)];
end
