function [meanY, E, pc] = calcpcStatic_Self(distMat, numTrials, numTargets, tempNum)

% This version chooses a single template per target for each iteration
numPerms = 100;

s = size(distMat); % [numTrains, numTrains, numTaus]

pc = zeros(prod(s(3:end)), s(1)/numTrials/numTargets, numTrials);
for permNum = 1:numPerms
    pc(:,:,permNum) = calcpcNew(distMat, numTrials, numTargets, tempNum);
end

E = std(pc,[],3);   % standard deviation
meanY = mean(pc,3);
meanY = reshape(meanY, [s(3:end), size(meanY, ndims(meanY)), 1]);
E = reshape(E, [s(3:end), size(E, ndims(E)), 1]);

%% the new shiny version
function pc = calcpcNew(distMat, numTrials, numTargets, tempNum)
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

inds1 = ceil((1:(numTrials-2)*numTargets)/(numTrials-2));    


% inds 2 = offset for distance matrix indexes(?)
inds2 = (0:numTrials:(numTargets - 1)*numTrials);

pc = zeros(numDistMats, 1);

for ii = 1:numDistMats
    tempMat = distMat(:,:,ii);    % for each tau value
    
    % reset(RandStream.getDefaultStream,0); % this can have negative effects on the code that calls it
    
    % choose random template for each target
    tempInds = floor(rand(1, numTargets)*numTrials);
    
    if all(tempInds == tempInds(1))
        tempInds(2) = mod(tempInds(2) + ceil(rand*(numTrials - 1)),numTrials); 
    end
    
    % top half of compInds for target 1 as trial, bottom half for target 2
    % as trial
    compInds = repmat(setxor([0:(numTrials-1)]',tempInds),numTargets,1);
    
    compInds = compInds(:) + 1; % shift compInds by 1 for indexing
    
    compShiftSingle = inds2;
    tempShiftSingle = circulant(inds2) + (tempNum - 1)*numTrials; % = [inds2, flipud(inds2)]
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

