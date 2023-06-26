%CALCPC Generate percent correct curves
%   [PC, E] = CALCPC(DISTMAT, NUMTRIALS, NUMTARGETS, TEMPNUM) uses the stacked
%   distance matrices in DISTMAT (a 3-dimensional array) to calculate mean
%   percent correct curves PC and the standard error E.  NUMTRIALS is the
%   number of trials recorded for each stimulus.  NUMTARGETS is the number of
%   choices the neuron has to classify a trial--it is equivalent to the
%   reciprocal of chance.  TEMPNUM is the index of the stimulus to use as the
%   template.
%
%   For example, if there are 6 total stimuli, 2 targets, and 10 trials of each
%   stimulus, one layer of the distance matrix will be 60 x 60.  If the middle
%   stimuli (the second and fifth blocks in the distance matrix) are to be used
%   as the template, TEMPNUM should be set to 2.
%
%   It is important to note that unlike the previous script, the distance
%   matrices do not need to be reordered.  Also, the percent correct values
%   are returned in the same order as the distance matrix, rather then the
%   template-first order as before.
%
%   It should also be noted that CALCPC is run after the distance matrices have
%   been created, which means it is neutral to the dimension over which the
%   different ditance matrices were made.  It could be tau, T, or anything
%   else.
%
%   [PC, E] = CALCPC(DISTMAT, NUMTRIALS, NUMTARGETS, TEMPNUM, NUMPERMS) is the
%   same as above, but specifies a number of permutations to perform.  Larger
%   numbers will lead to smaller values in E, but linearly increasing
%   computational time.

function [meanY, E, pc] = calcpcStatic(distMat, numTrials, numTargets, selfFlag)

% This version chooses a single template per target for each iteration and
% cycles through all possible pairs of templates
if ~selfFlag
    numIterations = numTrials ^ 2;
else
    tempMatrix = nchoosek(0:(numTrials-1),2);
    numIterations = size(tempMatrix,1);
end

s = size(distMat); % [numTrains, numTrains, numTaus]

for iterationNum = 1:numIterations
    if ~selfFlag
        [tempInds(1),tempInds(2)] = ind2sub([numTrials numTrials],iterationNum);
        tempInds = tempInds - 1;
    else
        tempInds = tempMatrix(iterationNum,:);
    end
    pc(iterationNum,:) = calcpcNew(distMat, numTrials, numTargets, tempInds, selfFlag);
end

E = std(pc);   % standard deviation
% maxY = mean(pc);
meanY = mean(pc);
% meanY = reshape(meanY, [s(3:end), size(meanY, ndims(meanY)), 1]);
% E = reshape(E, [s(3:end), size(E, ndims(E)), 1]);

%% the new shiny version
function pc = calcpcNew(distMat, numTrials, numTargets, tempInds, selfFlag)
if exist('calcpcHelper', 'file') == 3
	yFun = @calcpcHelper;
else
	yFun = @calcpcHelperSlow;
end

% tempNum will always be equal to 1 if you only have one unique response
% per target. if there were multiple responses for each target (for
% example, what if distMat contained the distances between a clean response
% and a masked response for each target?), tempNum will have to be changed

s = size(distMat);
n = s(1); % number of spike trains
numDistMats = prod(s(3:end)); % corresponds to number of tau values

% inds1 = which target response is being used for comparison between
% templates

if ~selfFlag   % comparing two different target responses = (numTargets * (numTrials - 1) matches)
    inds1 = ceil((1:(numTrials-1)*numTargets)/(numTrials-1));
else   % comparing target response to itself = (numTargets * (numTrials - 2) matches)
    inds1 = ceil((1:(numTrials-2)*numTargets)/(numTrials-2));
end

% inds 2 = offset for distance matrix indexes(?)
inds2 = (0:numTrials:(numTargets - 1)*numTrials);

pc = zeros(numDistMats, 1);

for ii = 1:numDistMats
    tempMat = distMat(:,:,ii);    % for each tau value
    
    compInds = [];
    
    % top half of compInds for target 1 as trial, bottom half for target 2
    % as trial
    if ~selfFlag
        for ti = 1:numTargets
            compInds(:,ti) = setxor(0:(numTrials-1),tempInds(ti));
        end
    else
        compInds = repmat(setxor([0:(numTrials-1)]',tempInds),numTargets,1);
    end
    
    compInds = compInds(:) + 1; % shift compInds by 1 for indexing
    
    compShiftSingle = inds2;
    tempShiftSingle = circulant(inds2); % = [inds2, flipud(inds2)]
    compShift = compShiftSingle(inds1).';
    tempShift = tempShiftSingle(inds1,:);
    
    a = compInds + compShift; % shift compInds for comparisons where responses to the 2nd target are the input
    
    for ti = 1:numTargets
        b = a + (tempInds(:,ti) + tempShift(:,ti))*n; % convert tempInds and tempShift to linear indices in tempMat
        distances(:,ti) = tempMat(b); % get indexes of the distance between the given template and the (ti)th target
    end
    
    % 1st column of distances will always be the correct target (1st 1000
    % rows have target 1 in column 1, other 1000 rows have target 2 in
    % column 1)
    
    % correct == 1 corresponds to where distances(n,1) > distances(n,2)
    
    correct = yFun(distances);
    pc(ii) = 100*mean(correct);
    % stdev(ii) = 100*std(correct, 1)/sqrt(numTargets);
end

function correct = calcpcHelperSlow(distances)
[minVal, order] = min(distances, [], 2);
correct = order(:,1) == 1; % these are only potentially correct (could be ties)
correct(correct) = rand(sum(correct), 1) < 1./sum(distances(correct,:) == repmat(minVal(correct), 1, size(distances,2)), 2); % now break any ties

