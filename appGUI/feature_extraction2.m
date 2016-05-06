function [ feature_mtx, label_mtx ] = feature_extraction2( preprocess_images , order )
%% STEP 4: Extracting Features From Each Pre-Processed Sample
label_mtx = [];
feature_mtx = [];
for ii = 1:length(preprocess_images)
    label_mtx = [label_mtx; ii*ones(size(preprocess_images{ii},3),1)];
    for jj = 1:size(preprocess_images{ii},3)
        [ xxx num_island] = bwlabel(preprocess_images{ii}(:,:,jj));
        feature_vector = [num_island];
        boundary_images = bwboundaries(preprocess_images{ii}(:,:,jj));
        boundary_cat_images = [];
        for kk = 1:length(boundary_images)
            boundary_cat_images = [boundary_cat_images; boundary_images{kk}];
        end
        rFSDs = fEfourier(boundary_cat_images, order,0,1);
        feature_vector = [feature_vector reshape(rFSDs,1,[])];
        feature_mtx = [feature_mtx;feature_vector];
    end
end

%%save('feature_mtx.mat','feature_mtx');
%%save('label_mtx.mat','label_mtx');
end