function [decision_2pin, decision_3pin] = KNN_sorted (feature_mtx_dir,label_mtx_dir, components_preprocess, pins);
%% Sort Then Classify
% Sort as 2pin or 3pin device
% First Sort Training Data
feature_mtx = importdata(feature_mtx_dir);
label_mtx = importdata(label_mtx_dir);
%preprocessed_images = importdata('preprocess_images.mat');
%[feature_mtx,label_mtx] = feature_extraction2(preprocessed_images, lol*5);

feature_2pin = []; label_2pin = [];
feature_3pin = []; label_3pin = [];
for ii = 1:length(label_mtx)
    xx = label_mtx(ii);
    if (xx == 1 || xx==2 || xx==3 ||  xx==10) %No such thing as inductor
        feature_2pin = [feature_2pin; feature_mtx(ii,:)];
        label_2pin = [label_2pin; label_mtx(ii)];
    elseif (xx == 6 || xx==7 || xx==8 || xx==9)
        feature_3pin = [feature_3pin; feature_mtx(ii,:)];
        label_3pin = [label_3pin; label_mtx(ii)];
    end
end

%% Second Sort Test Data
%components_preprocess = importdata('components_preprocess.mat');
for ii = 1:size(components_preprocess{1},3)
    %figure; imshow(components_preprocess{1}(:,:,ii));
    image = components_preprocess{1}(:,:,ii);
    yy = regionprops(image,'BoundingBox');
    cropped_image = imcrop(image, yy.BoundingBox);
    rescaled_image =imresize(cropped_image,[250 250]);
    thin_image = bwmorph(rescaled_image,'thin',Inf);
    components_preprocess{1}(:,:,ii) = thin_image;
    %figure;imshow(thin_image);
end

test_mtx = feature_extraction(components_preprocess);
test_2pin=[];
test_3pin=[];
for ii = 1:length(pins)
    if (pins(ii) == 2)
        test_2pin = [test_2pin; test_mtx(ii,:)];
    elseif (pins(ii) == 3)
        test_3pin = [test_3pin; test_mtx(ii,:)];
    end
end

%% Last, Classify
if (~isempty(test_2pin))
    decision_2pin = KNN_classifier(feature_2pin,label_2pin,test_2pin);
    decision_2pin = decision_2pin{1};
else
    decision_2pin = [];
end
if(~isempty(test_3pin))
    decision_3pin = KNN_classifier(feature_3pin,label_3pin,test_3pin);
    decision_3pin = decision_3pin{1};
else
    decision_3pin = [];
end

end