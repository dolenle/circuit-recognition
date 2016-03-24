function [ training_images,preprocess_images,feature_mtx,label_mtx ] = labeled_feature_mtx( directory, destination, guiThing )
%% Getting Training Images from Raw Scanned Images
tic
set(guiThing.textStatus, 'String', 'Loading training images...'); drawnow;
directory = [directory  '\'];
training_images = getting_training_images(directory);
%% Preprocessing Training Images
set(guiThing.textStatus, 'String', 'Preprocessing training images...'); drawnow;
preprocess_images = preprocess(training_images, guiThing); 
%% Feature Extraction From Preprocessed Images
set(guiThing.textStatus, 'String', 'Extracting Features...'); drawnow;
[feature_mtx label_mtx] = feature_extraction( preprocess_images);
%% Saving Matricies
set(guiThing.textStatus, 'String', 'Saving Data...'); drawnow;
save(fullfile(destination, 'training_images.mat'),'training_images')
save(fullfile(destination, 'preprocess_images.mat'),'preprocess_images');
save(fullfile(destination, 'feature_mtx.mat'),'feature_mtx');
save(fullfile(destination, 'label_mtx.mat'),'label_mtx');
set(guiThing.textStatus, 'String', ['Done. ' num2str(toc) 's elapsed.']); drawnow;
end
