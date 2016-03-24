%clear all; close all;
%% Importing Data
x_train = importdata('feature_mtx.mat');
y_train = importdata('label_mtx.mat');
x_test = importdata('feature_mtx.mat');
%% Dimensionality Reduction
confidence = 0.999;     %Correlation threshold
cov_mtx = cov(x_test);  %Covariance matrix of features
cor_mtx = corrcov(cov_mtx); %Correlation matrix of features
cor_mtx_lower = tril(cor_mtx-triu(cor_mtx));    %Take only the lower triangular elements of the correlation matrix
num_features = size(cor_mtx,1);                 %Number of features in original dataset
stop = 0;
redundant_features_test = zeros(1,num_features);
while stop == 0
    correlate = abs(cor_mtx_lower)>confidence;
    correlate = find(correlate);
    if isempty(correlate)
       stop = 1;
    else  
       row = mod(correlate(1)-1,num_features) + 1;
       column = (correlate(1)-row)/num_features + 1;
       cor_mtx_lower(row,:) = 0;
       cor_mtx_lower(:,row) = 0;
       redundant_features_test(row) = 1;
    end
end
cov_mtx = cov(x_train);   %Covariance matrix of features
cor_mtx = corrcov(cov_mtx); %Correlation matrix of features
cor_mtx_lower = tril(cor_mtx-triu(cor_mtx));  %Take only the lower triangular elements of the correlation matrix
num_features = size(cor_mtx,1);
stop = 0;
redundant_features_train = zeros(1,num_features);
while stop == 0
    correlate = abs(cor_mtx_lower)>confidence;
    correlate = find(correlate);
    if isempty(correlate)
       stop = 1;
    else  
       row = mod(correlate(1)-1,num_features) + 1;
       column = (correlate(1)-row)/num_features + 1;
       cor_mtx_lower(row,:) = 0;
       cor_mtx_lower(:,row) = 0;
       redundant_features_train(row) = 1;
    end
end
redundant_features = redundant_features_test .* redundant_features_train;   %Redundant features from training set AND test set

number_of_features_removed = sum(redundant_features)    %Number of extraneous features for given threshold 'confidence'
reduced_x_train = [];
reduced_x_test = [];

%Removing redundant columns from x_train and x_test 
for ii = 1:size(x_train,2)
    if redundant_features(ii) == 0
        reduced_x_train = [reduced_x_train x_train(:,ii)];
        reduced_x_test = [reduced_x_test x_test(:,ii)];
    end
end
%% Re-sorting Data into their Classes
% p is a number between 0 and 1 representing the percent of the training
% data used as a pseudo test_data.
% p was usually selected to be 0.33 such that there was a 2:1 ratio between
% training:psuedo test.
% When p == 0, the actual classification on 'test.txt' is performed using a
% pruned training data set.
p = 0.9;

yx_train = [y_train reduced_x_train];
yx_sort = sortrows(yx_train,1); %Groups all the categories from 1 to 12

%Retrieving index information to parse yx_sort
number_per_cat = zeros(12,1);
cat_index = zeros(12,1);
for ii=1:12
    number_per_cat(ii) = sum(ii==y_train);
    cat_index(ii)= sum(number_per_cat);
end

%Parsing the data into training and testing data in a cell corresponding to
%the actual label
train_data = cell(1,12);    % (1-p) of the data will be used for training
test_data = cell(1,12);     % (p)   of the data will be used for psuedo_testing
for ii = 1:12
    if ii == 1
        train_data{ii} = yx_sort(1:cat_index(ii)-floor(p*number_per_cat(ii)),2:end);
    else
        train_data{ii} = yx_sort(cat_index(ii-1)+1:cat_index(ii)-floor(p*number_per_cat(ii)),2:end); 
    end
    test_data{ii} = yx_sort(cat_index(ii)-floor(p*number_per_cat(ii))+1:cat_index(ii),2:end);
end

% In the case of classifying the actual data, the test_data is
% dimension reduced_x_test.
% The training set is a pruned, dimension reduced train_data
if p == 0;
    test_data{1} = reduced_x_test;
    train_data_pruned = importdata('train_data_pruned.mat');
    train_data = train_data_pruned;
end
%% Weighted K-Nearest Neighbors (with known labels)
tic
k = 13;     %k parameter for k-nearest neighbors; it is better for k to be odd and small in value

score = zeros(1,12);
decision = cell(1,12);

ii_end = 12;
if p == 0
    ii_end = 1;
end

for ii=1:ii_end % ii'th test category
    decision{ii} = zeros(1,size(test_data{ii},1));
    for jj=1:size(test_data{ii},1)  %jj'th row of test data in ii'th category
        distances = []; %reset distance vector 
        min_dist = zeros(k,2);
        for kk=1:12 %kk'th train category
            for mm=1:size(train_data{kk},1) %%mm'th row of train data in kk'th category
                distances = [distances ; [sum(abs(test_data{ii}(jj,:)-train_data{kk}(mm,:))) , kk]]; %Manhattan Distance (L1 norm)
                %distances = [distances ; [pdist([test_data{ii}(jj,:);train_data{kk}(mm,:)]) , kk]]; %Euclidean Distance (L2 norm)
            end
        end
        distances = sortrows(distances,1);  %Sort all distances
        min_dist = distances(1:k,:);        %Find k-nearest neightbors
        for nn=1:12
            %score(nn) = sum((min_dist(:,2) == nn)  %Uniform weighting
            %score(nn) = sum((min_dist(:,2) == nn)./(min_dist(:,1)+10e-12));   %1/distance weighting 
            %score(nn) = sum((min_dist(:,2) == nn).*(1+((min_dist(:,1)-min_dist(1,1))/(min_dist(1,1)-min_dist(k,1)))));  %Linear weighting normalized between 0 and 1
            score(nn) = sum((min_dist(:,2) == nn).*normpdf([min_dist(:,1)],0,min_dist(k,1)/3)); %Gaussian Weighting
        end
        [xxx decision{ii}(jj)] = max(score);    % Category bin with highest score wins
    end
end
%% Analysis and Pruning for p ~= 0
results = cell(1,12);
total_correct = 0;
total = 0;
test_data_pruned = cell(1,12);
train_data_pruned = cell(1,12);
for pp = 1:12
    results{pp} = sum(decision{pp} == pp)/length(decision{pp});
    test_data_pruned{pp} = test_data{pp}(decision{pp} == pp,:);
    train_data_pruned{pp} = [train_data{pp} ; test_data_pruned{pp}];
    total_correct = total_correct + sum(decision{pp} == pp);
    total = total + length(decision{pp});
end
toc

if p~=0
    results     % Percent Correct in each category
    percent_correct = total_correct/total   %Total percent correct
end



%% Outputting to .txt & .csv file
if p ==0
    fileID = fopen('Nguyen_Chowdhury_Output.txt','w');
    fprintf(fileID,'%d, %d\n', [(1:length(decision{1})) ; decision{1}]);
    csvwrite('Nguyen_Chowdhury_Output.csv',[(1:length(decision{1}))' , decision{1}']);
end

% Save 'train_data_pruned' variable as actual test data
if p ~= 0
    save('train_data_pruned.mat','train_data_pruned');
end