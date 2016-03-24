function [ feature_mtx, label_mtx ] = feature_extraction( preprocess_images )
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
%%save('feature_mtx.mat','feature_mtx');
%%save('label_mtx.mat','label_mtx');
end

