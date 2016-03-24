function [ training_images ] = getting_training_images( directory )
%% STEP 1: Getting Raw Scanned Images From Folders
% raw_images is a cell array of cell arrays containing the orginal
% scanned images
disp(directory);
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
        [columns, rows] = getGrid(14, 11, bw_image);
        columns = unique(columns');
        rows = unique(rows');
%         % Extract Row Indexes
%         rows = [];
%         for kk = 1:size(bw_image,1)
%             if (sum(bw_image(kk,:)) < (size(bw_image,2)*0.7))
%                 rows = [rows kk];
%             end
%         end
%         rows = rows(([250 (rows(2:end) - rows(1:end-1))] > 2));
%         % Extract Column Indexes
%         columns = [];
%         for kk = 1:size(bw_image,2)
%             if (sum(bw_image(:,kk)) < (size(bw_image,1)*0.75))
%                 columns = [columns kk];
%             end
%         end
%         columns = columns(([250 columns(2:end) - columns(1:end-1);] > 5));
%         % Hard Coded Row and Column Indexes for NON-CAPACITORS
%         if ii ~= 2
%             rows = [rows(1) rows(1)+265*[1:13]];
%             columns = [columns(1) columns(1)+263*[1:10]];
%         end
%         ericrow=rows
%         ericcol = columns
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
%save('training_images.mat','training_images');
end