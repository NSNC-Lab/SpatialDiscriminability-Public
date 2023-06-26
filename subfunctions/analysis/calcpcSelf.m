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

function [meanY, E, correctArray, pc] = calcpcSelf(distMat, numTrials, numTargets, tempNum, numPerms, ver)

if ~exist('ver', 'var') || isempty(ver)
	ver = 'new';
end
ver = lower(ver);

switch ver
    case 'new'
        if ~exist('numPerms', 'var') || isempty(numPerms)
            numPerms = 1e4;
        end
		s = size(distMat); % [numTrains, numTrains, numTaus]
% 		pc = zeros(prod(s(3:end)), s(1)/numTrials/numTargets, numTrials);        
% 		for trialNum = 1:numTrials
% 			inds = setxor(trialNum:numTrials:s(1), 1:s(1));
%             %C = setxor(A,B) for vectors A and B, returns the values that 
%             %are not in the intersection of A and B with no repetitions. 
%             %C will be sorted
% 			pc(:,:,trialNum) = calcpcNew(distMat(inds, inds, :), numTrials - 1, numTargets, tempNum, ceil(numPerms/numTrials));
% 		end
%         
%         % LOO approach for all possible pairs of removed trials
%         pc = zeros(prod(s(3:end)), s(1)/numTrials/numTargets, numTrials^2);        
% 		for trialNum = 1:numTrials^2
%             [r1,r2] = ind2sub([numTrials numTrials],trialNum);
% 			inds_t1 = setxor(r1:numTrials:s(1), 1:s(1));
%             inds_t2 = setxor(r2:numTrials:s(1), 1:s(1));

        pc = zeros(prod(s(3:end)), s(1)/numTrials/numTargets, numTrials);
        for trialNum = 1:numTrials
            inds_t1 = setxor(trialNum:numTrials:s(1), 1:s(1));
            inds_t2 = setxor(trialNum:numTrials:s(1), 1:s(1));

			pc(:,:,trialNum) = calcpcNew(distMat(inds_t1, inds_t2, :), numTrials - 1, numTargets, tempNum, ceil(numPerms/numTrials));
		end

		E = std(pc, [], 3)/sqrt(size(pc,3));   % standard error
		meanY = mean(pc, 3);
	case 'old'
        if ~exist('numPerms', 'var') || isempty(numPerms)
            numPerms = 100;
        end
        [meanY, E] = calcpcOld(distMat, numTrials, numTargets, tempNum, numPerms);
    otherwise
        error('VER must me ''old'' or ''new''')
end
meanY = reshape(meanY, [s(3:end), size(meanY, ndims(meanY)), 1]);
E = reshape(E, [s(3:end), size(E, ndims(E)), 1]);
correctArray = [];

%% the new shiny version
function pc = calcpcNew(distMat, numTrials, numTargets, tempNum, numPerms)
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
numVersions = n/numTrials/numTargets; % corresponds to numMaskers

inds1 = ceil((1:numPerms*numTargets)/numPerms);
inds2 = (0:numVersions*numTrials:(numTargets - 1)*numVersions*numTrials);

pc = zeros(numDistMats, numVersions);
distances = zeros(numPerms*numTargets, numTargets);
for ii = 1:numDistMats
    tempMat = distMat(:,:,ii);
    for versionNum = 1:numVersions   % numVersions = 1 for our data
        % reset(RandStream.getDefaultStream,0); % this can have negative effects on the code that calls it
        
        compInds = floor(rand(numPerms*numTargets, 1)*numTrials);   % indexes of spike trains to match to templates
        tempInds = floor(rand(numPerms*numTargets, numTargets)*numTrials);  % indexes of templates for each comparison
        
        % shift tempInds by at least 1 index to prevent comparisons where
        % compInds = tempInds(:,1), which would always be correct
        tempInds(:,1) = mod(compInds + ceil(rand(numPerms*numTargets, 1)*(numTrials - 1)), numTrials);
        tempInds(:,2) = mod(compInds + ceil(rand(numPerms*numTargets, 1)*(numTrials - 1)), numTrials);
        
        compInds = compInds + 1; % shift compInds by 1 for indexing
        
        compShiftSingle = inds2 + (versionNum - 1)*numTrials; % =inds2 for our purposes, since numVersions = 1
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
        pc(ii,versionNum) = 100*mean(correct);
        % stdev(ii,versionNum) = 100*std(correct, 1)/sqrt(numPerms*numTargets);
    end
end

function correct = calcpcHelperSlow(distances)
[minVal, order] = min(distances, [], 2);
correct = order(:,1) == 1; % these are only potentially correct (could be ties)
correct(correct) = rand(sum(correct), 1) < 1./sum(distances(correct,:) == repmat(minVal(correct), 1, size(distances,2)), 2); % now break any ties

%% the old slow version
function calcpcOld(distMat, numTrials, numTargets, tempNum, numPerms)
nTests = size(distMat, 3);

numMaskers = size(distMat, 1)/numTargets/numTrials;
totPerm = numPerms*numTrials;

tempOrder = [tempNum, 1:tempNum - 1, tempNum + 1:numMaskers];
order = zeros(1, numMaskers*numTargets);
% order = zeros(1, numMaskers*2);
if numMaskers == 1
	order = [tempNum, 1:tempNum - 1, tempNum + 1:numTargets];
else
	for i = 1:numTargets
		order(i:numTargets:numMaskers*numTargets) =...
			tempOrder + numMaskers*(i-1);
	end
end
% order(1:2:numMaskers*2 - 1) = tempOrder;
% order(2:2:numMaskers*2) = tempOrder + numMaskers;
distMat = reorder_distmat(distMat, order);

% discardTrials = ud.spikes.record{1}.discardTrials;
% keepTrials = ones(ud.spikes.record{1}.numSpikeTrainsPerSong*ud.songs.header.numRecords,1);
% keepTrials(discardTrials) = 0;
  
% generate template permutations
templates = zeros(numTargets, totPerm);
for j = 1:numPerms
	for i = 1:numTargets
		colStart = ((j - 1)*numTrials) + 1;
        tpl = randperm(numTrials) +  numTrials*(i - 1);
% 		templates(i,colStart:colStart + numTrials - 1) = randperm(numTrials) +  numTrials*(i - 1);

        %Check for discarded trials, replace with valid trials
%         badIdx = find(~keepTrials(tpl));
%         while (~isempty(badIdx))
%             newTpl = randint(length(badIdx),1,[1,...
%                 ud.spikes.record{1}.numSpikeTrainsPerSong]);
%             tpl(badIdx)=newTpl;
%             badIdx = find(~keepTrials(tpl));         
%         end
        templates(i,colStart:colStart + numTrials - 1) = tpl;
	end
end

correctArray = zeros(nTests, totPerm, numMaskers);

for jTpl = 1:totPerm
	for iT = 1:nTests
		for iM = 1:numMaskers
			startidx = 1 + (iM - 1)*numTrials*numTargets;
			stopidx = startidx + (numTrials*numTargets) - 1;
			
			% update the confusion matrix
			confusionMatrix = confmat(distMat(startidx:stopidx,1:numTargets*numTrials,iT), templates(:,jTpl), numTargets,numTrials);

			% get the classification percent correct
			correctArray(iT,jTpl,iM) = 100*sum(diag(confusionMatrix))/((numTargets*numTrials) - numTargets);
		end
	end
end
for iM = 1:numMaskers
	meanY(:,order(numTargets*(iM - 1)+1)) = mean(correctArray(:,:,iM),2);
	E(:,order(numTargets*(iM - 1)+1)) = std(correctArray(:,:,iM),1,2);
end


%------------------------------------------------------------------------------%

%confmat.m
%compute the confusion Matrix given the distance matrix and a set of templates
% CM = CONFMAT(distMatrix, templates , numSongs, totSpikes, numSpikeTrains)
% distMatrix - distance matrix
% templates - sorted array of indices in distMatrix which are to be used as templates
% numSongs - number of songs
% numSpikeTrains - number of spike trains per song
% CM - confusion matrix (numSongs x numSongs)
% CM = CONFMAT(distMatrix, templates , numSongs, totSpikes, numSpikeTrains,excludeFlag)
% CM = CONFMAT(distMatrix, templates , numSongs, totSpikes, numSpikeTrains,excludeFlag,keepTrials)
% excludeFlag - include (0) or exclude(1,default) templates from
% being classified
% keepTrials - 
%Created 6/13/2003
%Changes: 6/16/2004 - excluded templates from classification
%                   - ties are now resolved by randomly assigning
%                   the winning template  
%         11/9/2005 - added an exclude_flag to either include or
%         exclude the templates from being classified
%         6/19/2006 - added code to deal with discarded trials
%Author: Rajiv Narayan

function confMatrix=confmat(distMatrix,templates,numSongs,numTrials,varargin)

%we can choose to exclude / not exclude the templates from
%classification.

exclude_flag=1; %default is to exclude the templates
keepTrials = ones(numSongs*numTrials,1); %default is to keep all trials

if (nargin>5)
	keepTrials = varargin{2};
end

if (nargin>4)
	exclude_flag=varargin{1}; %set the exclude_flag specified by user
else

end

if (nargin>=4)
	totSpikes = numSongs * numTrials;

	%Allocate space for the confusion Matrix
	confMatrix=zeros(numSongs);
	for  i=1:totSpikes

		if (keepTrials(i)) %valid trial
			tplcount = fix((i - 1) / numTrials) + 1;
			if (exclude_flag)
				%exclude templates from classification
				if (~isequal(templates(tplcount),i))
					%find the minimum dist
% 					[minDist, assignedSong] = min(distMatrix(i,templates(randperm(tplcount))));
					minDist=min(distMatrix(i,templates));

					%locate indices of ties if any
					tieList=find(distMatrix(i,templates)==minDist);
					numties=length(tieList);
					%we randomly assign the song to one of the ties
					assignedSong=tieList(round(rand*(numties-1))+1);
					songNum =  fix((i - 1) / numTrials) + 1;
					confMatrix(songNum,assignedSong) = confMatrix(songNum,assignedSong) + 1;
				end
			else
				%find the minimum dist
				minDist=min(distMatrix(i,templates));

				%locate indices of ties if any
				tieList=find(distMatrix(i,templates)==minDist);
				numties=length(tieList);
				%we randomly assign the song to one of the ties
				assignedSong=tieList(round(rand*(numties-1))+1);
				songNum =  fix((i - 1) / numTrials) + 1;
				confMatrix(songNum,assignedSong) = confMatrix(songNum,assignedSong) + 1;
			end

		end

	end
else
	errordlg('Insufficient Arguments',mfilename,'replace');
end
