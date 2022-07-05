function [confMat meanAcc] = crossValClassify(data_rvc,labels_rc,k,rand)

%-------------------------------------------
%
% crossValClassify(data_rvc,labels_rc,k,rand)
%
% do a cross-validated classification analysis 
% on a set of responses from multiple voxels / neurons
% to multiple conditions across multiple runs
%
% required input:
% data_rvc -- data matrix, runs by signals by conditions
% labels_rc -- label matrix, runs by conditions, must be integers
%
% optional inputs:
% k -- number of cross-validation folds (default = 10)
% rand -- permute the labels? (default = 0)
%
% freeman, 04-05-2012
%-------------------------------------------

if ~exist('rand','var') || isempty(rand)
    rand = 0;
end

if ~exist('k','var') || isempty(k)
    k = 10;
end

% get the dimensions
R = size(data_rvc,1);
V = size(data_rvc,2);
C = size(data_rvc,3);

% get the cross validation folds
cv = cvpartition(1:R,'Kfold',k);

% create an empty confusion matrix
confMat = zeros(C,C,cv.NumTestSets);
for icv=1:cv.NumTestSets
    % get training and testing indices
    trIdx = cv.training(icv);
    teIdx = cv.test(icv);
    % store training and testing data
    train_data = reshape(data_rvc(trIdx,:,:),[cv.TrainSize(icv)*C V]);
    test_data = reshape(data_rvc(teIdx,:,:),[cv.TestSize(icv)*C V]);
    % store training and testing labels
    train_labels = reshape(labels_rc(trIdx,:,:),[cv.TrainSize(icv)*C 1]);
    test_labels = reshape(labels_rc(teIdx,:,:),[cv.TestSize(icv)*C 1]);
    % randomize training labels if we're doing a permutation test
    if rand
        train_labels = shuffle(train_labels);
    end
    % estimate the labels on testing data
    est_labels = classify(test_data,train_data,train_labels,'diagLinear');
    % fill in the confusion matrix
    for iC_1=1:C
        for iC_2=1:C
            confMat(iC_1,iC_2,icv) = sum(test_labels==iC_1 & est_labels==iC_2);
        end
    end
end

% normalize the confusion matrix
confMat = sum(confMat,3)/sum(cv.TestSize);
% get average accuracy across conditions
meanAcc = mean(diag(confMat));