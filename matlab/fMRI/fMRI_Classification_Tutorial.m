%% MRI Classification Tutorial ------------
% 
% by j freeman for CSHL, 2012
%
addpath('dependencies_fMRI_Classification_Tutorial');
%%-----------------------------------------
%%
% fMRI data are inherently multivariate. We typically measure the BOLD
% response to several different conditions (e.g. orientations, objects) in
% a large number of voxels. We want to establish a statistical dependency
% between the conditions and the voxel responses. 
%
% Traditional analyses ("mass-univariate") separately estimate the response
% of each voxel to each condition using a general linear model (see
% ER_fMRI_Tutorial.m). fMRI responses are noisy, so the response of any
% individual voxel will rarely distinguish robustly between even a pair of
% conditions (let alone several conditions). One way to combat noise is to
% average the responses of many voxels together. If a large continugous
% group of voxels all differentially respond to one condition over another,
% their mean response may reliably differ between the two conditions. This
% is the basis for "ROI" analyses, as well as mass-univariate analyses that
% incoroproate spatial contiguity.
%
% But if the preferences are distributed, with each voxel showing a
% slightly different bias for one condition or another, we need statistical
% tests that pool all of these signals together to establish a dependency
% between the condition and the *pattern* of responses. This is the basis
% for "multivariate analyses".
%
% For the purposes of this tutorial, we'll assume that we have already
% obtained a response amplitude to each condition from each voxel by
% performing standard event-related analyses. And we'll additionally assume
% that for each voxel and condition, we've obtained responses across
% multiple independent runs. We'll generate fake data with this form, show
% why univariate anlayses are ineffecitve, and show how to run a simple
% classification analyses and do some statistics on the output. Note that
% all of the techniques described below could be applied nearly identically
% to populations of neural responses -- see Rust and DiCarlo, Neuron, 2012
% for a useful review.
%%
clear all
clf

% set some constants for the data:
V = 50; % number of voxels (default = 50)
R = 10; % number of runs (default = 10)
C = 5; % number of conditions (default = 5)
sigma = 0.5; % noise of fMRI response amplitude estimates (default = 0.5)

% set some parameters for the analysis
k = 10; % number of folds for cross-validation
nBoot = 10; % number of bootstraps (anything above 100 will be too slow!)
nPerm = 100; % number of permutations for randomization test

%%
% Each condition evokes the same mean response, with a fixed but
% essentially arbitrary pattern of responses across the voxels
mnVal = 1; % mean response
pattern_nc = rand(V,C)+(mnVal-0.5); % voxel patterns

% The responses in each run are given by the fixed voxel-by-condition
% pattern plus additive gaussian noise The complete data matrix is RxVxC --
% runs by voxels by conditions
data_rvc  = zeros(R,V,C);
for iR=1:R
    data_rvc(iR,:,:) = pattern_nc + randn(V,C)*sigma; % add noise
end

%%
% Look at the voxel-by-run patterns for the first two conditions. The
% responses on each run should look like a noisy version of the common
% pattern. Try turning the noise (sigma) down to see more clearly that the
% two patterns are distinct.
subplot(3,4,1);
imagesc(data_rvc(:,:,1),[min(data_rvc(:)) max(data_rvc(:))]);
xlabel('Voxel');
ylabel('Run');
title('Condition 1');
subplot(3,4,2);
imagesc(squeeze(data_rvc(:,:,2)),[min(data_rvc(:)) max(data_rvc(:))]);
xlabel('Voxel');
ylabel('Run');
title('Condition 2');
colormap(gray);

% The mean response (across voxels and runs) does not distinguish well
% between the different conditions.
subplot(3,4,3);
mn = [squeeze(mean(mean(data_rvc,1),2))]; % mean response to each condition
err = [squeeze(std(mean(data_rvc,2)))]; % standard deviation across runs
bar(mn); hold on;
h = errorbar([1:C],mn,2*err,2*err,'k'); hold off;
set(h,'LineStyle','none');
ylim([0 max(mn)*1.2]);
xlabel('Condition');
ylabel('Mean response');
box off;
xlim([0 C+1]);

%%
% For a few individual voxels, look at the distribution of responses to the
% first two conditions. For some voxels the distributions are nearly
% identical, for others they are more separate, but in each case there is
% substantial overlap.
b_resp = 10; % bin width for response histograms
for iVox=1:3
    mnVal = min(vector(data_rvc(:,iVox,:)));
    mxVal = max(vector(data_rvc(:,iVox,:)));
    bins = linspace(mnVal,mxVal,b_resp);
    h1 = hist(data_rvc(:,iVox,1),bins); % histogram for condition 1
    h2 = hist(data_rvc(:,iVox,2),bins); % histogram for condition 2
    subplot(3,4,4+iVox);
    plot(bins,h1/sum(h1*1/b_resp),'b--'); hold on; 
    plot(bins,h2/sum(h2*1/b_resp),'r--'); hold off;
    xlabel(sprintf('Voxel %g response',iVox));
    ylabel('Probability');
    box off;
end
%%
% One way to quantify the distinctiveness of these distributions is to ask
% how well we could use the response distributions estimated on one set
% of data to predict condition identity from novel responses. To make that
% precise, we'll fit a simple model to each response distributions -- a
% gaussian -- using data from all but the last run. In 1-D, we just need
% the mean and standard deviation.
b_fit = 100; % bin width for gaussian evaluation
for iV=1:3
    mnVal = min(vector(data_rvc(:,iV,:)));
    mxVal = max(vector(data_rvc(:,iV,:)));
    x_vals = linspace(mnVal,mxVal,b_fit);
    mn_1(iV) = mean(data_rvc(1:R-1,iV,1)); % get the condition 1 mean
    std_dev_1(iV) = std(data_rvc(1:R-1,iV,1)); % get the condition 1 std
    h1_fit = normpdf(x_vals,mn_1(iV),std_dev_1(iV)); % evaluate the gaussian
    mn_2(iV) = mean(data_rvc(1:R-1,iV,2)); % same for condition 2...
    std_dev_2(iV) = std(data_rvc(1:R-1,iV,2));
    h2_fit = normpdf(x_vals,mn_2(iV),std_dev_2(iV));
    subplot(3,4,4+iV);
    hold on;
    plot(x_vals,h1_fit/sum(h1_fit*1/b_fit),'b');
    plot(x_vals,h2_fit/sum(h2_fit*1/b_fit),'r'); hold off
    box off;
end

%%
% To estimate our prediction accuracy, we take a novel response that was
% not used to fit the gaussians, and ask whether it's more likely under the
% first or second gaussian
fprintf('\r1-D classification... \r');
for iV=1:3
    % for each voxel, compute the probability of a novel response from
    % condition 1 under the two different gaussians
    prob_1 = normpdf(data_rvc(R,iV,1),mn_1(iV),std_dev_1(iV));
    prob_2 = normpdf(data_rvc(R,iV,1),mn_2(iV),std_dev_2(iV));
    % compare the probabilities (it's usually easier to compare the logs
    % because they are more numerically stable when probabilities become
    % small)
    if log(prob_1) > log(prob_2)
        fprintf('Voxel %g says condition 1\r',iV);  
    else
        fprintf('Voxel %g says condition 2\r',iV); 
    end
end
%%
% One voxel at a time wasn't very useful. Let's look at the *joint*
% response of two voxels. Fit a 2-D gaussian to each condition's response
% distribution. As before, we'll predict the condition of a novel *pair* of
% voxel responses by computing probability under the gaussians.
fprintf('\r2-D classification... \r');
voxPairs = [1 2; 1 3; 2 3];
for iVp=1:3
    subplot(3,4,8+iVp);
    vox1 = voxPairs(iVp,1);
    vox2 = voxPairs(iVp,2);
    plot(data_rvc(:,vox1,1),data_rvc(:,vox2,1),'.'); 
    hold on; 
    plot(data_rvc(:,vox1,2),data_rvc(:,vox2,2),'.r');  
    hold off
    mn_1 = mean(data_rvc(:,[vox1 vox2],1)); % condition 1 mean 
    cov_1 = cov(data_rvc(:,[vox1 vox2],1)); % condition 1 cov matrix
    mn_2 = mean(data_rvc(:,[vox1 vox2],2)); % condition 2 mean
    cov_2 = cov(data_rvc(:,[vox1 vox2],2)); % condition 2 cov matrix
    drawEllipse(mn_1,cov_1,'b'); % plot the 1 SD contours of the gaussians
    drawEllipse(mn_2,cov_2,'r');
    xlabel(sprintf('Voxel %g response',voxPairs(iVp,1)));
    ylabel(sprintf('Voxel %g response',voxPairs(iVp,2)));
    box off;
    % compare the probability of a novel response from condition 1 under
    % the two different gaussians
    prob_1 = mvnpdf(data_rvc(R,[vox1 vox2],1),mn_1,cov_1);
    prob_2 = mvnpdf(data_rvc(R,[vox1 vox2],1),mn_2,cov_2);
    if log(prob_1) > log(prob_2)
        fprintf('Voxel pair %g says condition 1\r',iVp);
    else
        fprintf('Voxel pair %g says condition 2\r',iVp);
    end
end
% Try running the tutorial up to this point multiple times. In many cases,
% the 2-D distributions should appear more reliably distinct than they did
% in 1-D. There are several possible behaviors here. If the biases for two
% individual voxels are small but present, they will be more
% distinguishable when looking at them jointly. If the bias is only present
% in one voxel or another, we will get the right answer by looking at both.
% We don't need to pick the *right* one ahead of time.
%
%%
% We actually have V voxels. Rather than just look at each one, or each
% pair, we would like to look at all at once. For our full data set, the
% set of responses (across runs) to a particular condition is a cloud of
% points in a V-dimensional space, where each dimension is a voxel.
% Different conditions yield different clouds. We cannot visualize this
% space, but we can use the same procedure as above to carve up the space
% and make predictions about novel responses -- fit a V-Dimensional
% gaussian to each condition, and compute the probability of each novel
% response under the different Gaussians. That's exactly what the matlab
% function "classify.m" does! Let's use it.
%
% It's important that we fit our gaussians ("train") using one set of
% responses, and make predictions ("test") using novel data. This is called
% cross-validation.

% Create a CxR matrix of condition labels for each condition and run
labels_rc = repmat([1:C],[R 1]);

% Use the very helpful built-in function "cvpartition.m" to split the data
% into multiple training / testing sets. We'll use k-fold cross-validation,
% we means we divide our runs into k distinct subsets. We train on k-1 and
% test on the held out set.
cv = cvpartition(1:R,'Kfold',k);

% We just made an object (cv) that contains vectors of indices for each
% training / testing set. Type cv.traning(1) and cv.testing(1) to see the
% run indices for the first training / testing set.

%%
% Create an empty matrix to store the classification output
confMat = zeros(C,C,cv.NumTestSets);

% Now loop over the k training / testing sets
for icv=1:cv.NumTestSets
    % get the indices for the training runs
    trIdx = cv.training(icv);
    % get the indices for the testing runs
    teIdx = cv.test(icv);
    % these are the responses from the training and testing runs
    train_data = reshape(data_rvc(trIdx,:,:),[cv.TrainSize(icv)*C V]);
    test_data = reshape(data_rvc(teIdx,:,:),[cv.TestSize(icv)*C V]);
    % and these are the labels
    train_labels = reshape(labels_rc(trIdx,:),[cv.TrainSize(icv)*C 1]);
    test_labels = reshape(labels_rc(teIdx,:),[cv.TestSize(icv)*C 1]);
    % the "classify.m" function takes as input a matrix of test data, a
    % matrix of training data, and a vector of labels. With the option
    % 'diagLinear', it uses the training data to fit the distribution of
    % responses in each condition with a multivariate gaussian, and then
    % evaluates the relative probability of each of the responses in the
    % test data in order to assign a label to them -- just as we did before
    % in 1-D and 2-D. look inside to see more details...
    est_labels = classify(test_data,train_data,train_labels,'diagLinear');
    % store the results in a "confusion matrix" that collects the fraction
    % of time each "true" condition was classified as one of the possible
    % estimated conditions
    for iC_1=1:C
        for iC_2=1:C
            confMat(iC_1,iC_2,icv) = sum(test_labels==iC_1 & est_labels==iC_2);
        end
    end
end

% average the confusion matrix across all the k cross-validation splits
confMat = sum(confMat,3)/sum(cv.TestSize);
% the "accuracy", or percent correct, is the frequency with which the true
% response was correctly categorized -- the diagonal of the confusion
% matrix
acc = mean(diag(confMat));

% plot the confusion matrix from the classification. if we are doing well,
% the bright squares should be along the main diagonal.
subplot(3,4,4);
imagesc(confMat,[0 1]);
title(sprintf('Accuracy is %g',acc));
xlabel('True class');
ylabel('Estimated class');

% IMPORTANT NOTE: With the 'diagLinear' option, the classify function fits
% an isotropic gaussian to each condition's response distribution -- an
% isotropic gaussian can have a different variance in each direction, but
% the covariances (correlations) are assumed to be 0. In general, we expect
% voxel responses to be correlated. But we rarely have enough data to
% reliably estimate the correlations. A classifier that ignores the
% correlations reflects a lower bound on performance (it will always do
% worse than a classifier that knew about the correlations).

%%
% We can use simple resampling and randomization procedures to obtain
% statistical information about our classifier.

% How does the classifier performance vary as a function of the number of
% voxels included in the analysis? Loop over all possible number of voxels,
% and for each, grab several different random combniations of voxels and do
% the classification. We'll use a wrapper "crossValClassify.m" that
% collects some of the above into a single function. 
% WARNING: This part might be slow
acc_boot = zeros(V,nBoot);
for iBoot=1:nBoot
    for iV=1:V
        % get a random sample of voxels
        rndInds = randperm(V); 
        rndInds = rndInds(1:iV);
        % do the classification (we just need the accuracy)
        [~, acc_boot(iV,iBoot)] = crossValClassify(data_rvc(:,rndInds,:),labels_rc,k);
    end
end
subplot(3,4,8);
plot([1:V],prctile(acc_boot,50,2),'k','LineWidth',4); hold on;
plot([1:V],prctile(acc_boot,[5 95],2),'k'); hold off;
drawHorzLine(1/C,'k',':');
ylim([0 1.1]);
xlim([0 V]);
xlabel('Number of voxels included');
ylabel('Accuracy');
box off;

% Often, we will be comparing classificaiton across multiple brain areas
% that differ in their voxel counts (or the number of neurons, for
% physiology data sets). The above clearly shows that we should compare
% classification accuracy only after matching for the number of voxels!
%%
% Is our classification reliably above chance? Our null hypothesis is that
% there is no relationship between the identity of each condition and the
% responses on that condition. Can we reject the null? To simulate the
% distribution of classification under the null hypothesis, we'll classify
% our data several times, each time randomly permuting the condition labels
% on the training data.
acc_permuted = zeros(1,nPerm);
for iPerm=1:nPerm
    [~, acc_permuted(iPerm)] = crossValClassify(data_rvc,labels_rc,k,1);
end
% compute a p value as the fraction of the null distribution greater than
% the accuracy we observed on the real (unpermuted) data
pVal = sum(acc<acc_permuted)/nPerm;
subplot(3,4,12);
bins = linspace(0,1,25);
h = hist(acc_permuted,bins);
bar(bins,h/sum(h)); shading flat;
drawVertLine(acc);
drawVertLine(1/C,'k',':');
xlim([-0.1 1.1]);
box off;
xlabel('Accuracy');
ylabel('Frequency');
title(sprintf('P value is %g',pVal))

%%
% Now go back and run the script with different parameters!

% If you jack up the noise (to 1) or the number of conditions (to 20), or
% both, you should find that classification accuracy starts to take a hit.
% But you can recover it somewhat by either increasing the number of runs
% (try 40) or increasing the number of voxels (try 300).