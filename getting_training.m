clc;clear all;close all; format compact; tic;
%% STEP 1: Getting Raw Scanned Images From Folders
% raw_images is a cell array of cell arrays containing the orginal
% scanned images
directory = 'C:\Users\Dolen\Documents\MATLAB\SrProject\training_data\';
training_data_dir = dir(directory);
subfolders = training_data_dir(3:end);
raw_images = cell(1,length(subfolders));
for ii = 1:length(subfolders)
    subfolder_dir = dir([directory subfolders(ii).name]);
    subfolder_dir = subfolder_dir(3:end);
    raw_images{ii} = cell(0,0);
    for jj = 1:length(subfolder_dir)
        raw_images{ii} = [raw_images{ii} imread(subfolder_dir(jj).name)];
    end
end
%% STEP 2: Extracting Individual Samples From Raw Scanned Images
width = 20;
training_images = cell(1,length(raw_images));
for ii = 1:length(raw_images)
    training_images{ii} = cell(length(raw_images{ii}),130);
    for jj = 1:length(raw_images{ii})
        % Retrieve One Binary image from raw_images
        bw_image = im2bw(raw_images{ii}{jj},0.95);
        % Extract Row Indexes
        rows = [];
        for kk = 1:size(bw_image,1)
            if (sum(bw_image(kk,:)) < (size(bw_image,2)*0.7))
                rows = [rows kk];
            end
        end
        rows = rows(([250 (rows(2:end) - rows(1:end-1))] > 2));
        % Extract Column Indexes
        columns = [];
        for kk = 1:size(bw_image,2)
            if (sum(bw_image(:,kk)) < (size(bw_image,1)*0.75))
                columns = [columns kk];
            end
        end
        columns = columns(([250 columns(2:end) - columns(1:end-1);] > 5));
        % Hard Coded Row and Column Indexes for NON-CAPACITORS
        if ii ~= 2
            rows = [rows(1) rows(1)+265*[1:13]];
            columns = [columns(1) columns(1)+263*[1:10]];
        end
%         % Deleting Rows and Columns
%         for kk = -width:width
%             bw_image(rows+kk,:) = 1;
%             bw_image(:,columns+kk) = 1;
%         end
%         figure; imshow(bw_image);

        % Cutting Out Individual Samples from BW-image
        for mm = 1:13
            for nn = 1:10
                a = rows(mm)+width;
                b = columns(nn)+width;
                c = rows(mm+1)-width;
                d = columns(nn+1)-width;
                training_images{ii}{jj,10*(mm-1)+nn} = ~bw_image(a:c,b:d);
            end
        end
        
    end
end
save('training_images.mat','training_images');
%% TEST 2.1: Check If Individual Samples Were Extracted Correctly
% close all;
% figure
% subplot(10,13,1)
% a = 1;
% b = 0;
% for ii = 1:130
%     subplot(13,10,ii);
%     imshow(training_images{a}{ii+b});
% end
%% STEP 3: Pre-Processing Each Sample
preprocess_images = cell(1,length(training_images));
for jj = 1:length(training_images)
    for kk = 1:(size(training_images{jj},1)*size(training_images{jj},2))
        % Retrieve sample from training_images
        sample = training_images{jj}{kk};
        
        %Denoising
        sample_dilation_denoised = imdilate(imerode(sample,strel('disk',1)),strel('disk',1));
        [sample_label num_islands] = bwlabel(sample_dilation_denoised);

        for aa = 1:num_islands
            if( sum(sum(sample_label == aa)) <= 125 )
                sample_label(sample_label == aa) = 0;
            end
        end
        sample_denoised = im2bw(sample_label);
    
        %Cutting & Rescaling
        ii=1;
        while((sum(sample_denoised(ii,:))) == 0)
            a = ii;
            ii = ii+1;
        end
        ii=0;
        while((sum(sample_denoised(end-ii,:)))==0)
            b = size(sample_denoised,1)-ii;
            ii = ii+1;
        end
        ii=1;
        while((sum(sample_denoised(:,ii))) == 0)
            c = ii;
            ii = ii+1;
        end
        ii=0;
        while((sum(sample_denoised(:,end-ii)))==0)
            d = size(sample_denoised,2)-ii;
            ii = ii+1;
        end

        sample_rescaled = imresize(sample_denoised(a+1:b-1,c+1:d-1),[250 250]);
        
        %Skeltonize
        %sample_skel = bwmorph(sample_rescaled,'skel',Inf);
        sample_skel = bwmorph(imdilate(sample_rescaled,strel('disk',3)),'thin',Inf);

        preprocess_images{jj} = cat(3,preprocess_images{jj}, sample_skel);

    end
end
save('preprocess_images.mat','preprocess_images');
%% TEST 3.1: Check If Individual Samples Were Pre-Processed Correctly
% close all;
% for ee = 1:10
%     figure; imshow(preprocess_images{10}(:,:,ee));
% end
%% TEST 3.2: Find the Average Number of Islands Per Class
% total = zeros(1,10);
% for ii = 1:10
%     for jj = 1:size(preprocess_images{ii},3)
%         [xxx yyy] = bwlabel(preprocess_images{ii}(:,:,jj));
%         total(ii) = total(ii) + yyy;
%     end
%     total_avg(ii) = total(ii)/size(preprocess_images{ii},3);
% end
% total_avg
%% STEP 4: Extracting Features From Each Pre-Processed Sample
% LIST OF PROPERTIES:
% 'Area','MajorAxisLength','MinorAxisLength','Centroid', ...
% 'Orientation','Eccentricity','ConvexArea','FilledArea', ...
% 'EquivDiameter','Extent','Solidity'
properties = {'Area','MajorAxisLength','MinorAxisLength' };

divisions = 5;
sub_index = floor(linspace(1,size(preprocess_images{1},1),divisions+1));
middle_sub_index = sub_index(2:end-1);
sub_index = reshape([sub_index(1) reshape([middle_sub_index;middle_sub_index+1],1,[]) sub_index(end)],2,[])';
[i_row j_col] = meshgrid(1:divisions);

label_mtx = [];
feature_mtx = [];
for ii = 1:length(preprocess_images)
    label_mtx = [label_mtx ii*ones(1,size(preprocess_images{ii},3))];
    for jj = 1:size(preprocess_images{ii},3)
        [ xxx num_island ] = bwlabel(preprocess_images{ii}(:,:,jj));
        feature_vector = [num_island];       
        for pp = 1:divisions.^2
            sub_preprocess_image = preprocess_images{ii}(sub_index(i_row(pp),1):sub_index(i_row(pp),2), sub_index(j_col(pp),1):sub_index(j_col(pp),2),jj);
            for kk = 1:length(properties)
                reg_prop = regionprops(sub_preprocess_image,properties{kk});
                property = struct2array(reg_prop);
                if size(property) == 0
                    property = 0;
                end
                feature_vector = [feature_vector property];
            end
        end
        feature_mtx = [feature_mtx;feature_vector];
    end
end
label_mtx = label_mtx';
save('feature_mtx.mat','feature_mtx');
save('label_mtx.mat','label_mtx');
%%
 toc;