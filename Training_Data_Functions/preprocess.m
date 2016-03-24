function [ preprocess_images ] = preprocess( training_images, guiHandle )
%% STEP 3: Pre-Processing Each Sample
set(guiHandle.textStatus, 'String', 'Preprocessing training images...'); drawnow;
totalSamples = 0;
counter = 0;
for jj = 1:length(training_images)
    totalSamples = totalSamples+size(training_images{jj},1)*size(training_images{jj},2);
end
totalSamples = num2str(totalSamples);
preprocess_images = cell(1,length(training_images));
for jj = 1:length(training_images)
    for kk = 1:(size(training_images{jj},1)*size(training_images{jj},2))
        counter = counter+1;
        set(guiHandle.textStatus, 'String', ['Preprocessing sample ',num2str(counter),' of ',totalSamples]); drawnow;
        
        % Retrieve sample from training_images
        sample = training_images{jj}{kk};
        axes(guiHandle.axes1);
        imshow(sample);
        
        %Denoising
        sample_dilation_denoised = imdilate(imerode(sample,strel('disk',1)),strel('disk',1));
        [sample_label num_islands] = bwlabel(sample_dilation_denoised);

        for aa = 1:num_islands
            if( sum(sum(sample_label == aa)) <= 125 )
                sample_label(sample_label == aa) = 0;
            end
        end
        sample_denoised = im2bw(sample_label);
    
%         box = regionprops(sample_denoised, 'boundingbox');
%         sample_rescaled = imresize(imcrop(sample_denoised, box(1).BoundingBox), [250 250]);
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
end

