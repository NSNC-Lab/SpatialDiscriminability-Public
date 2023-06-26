
function [clusterIDs, unitQuality, contaminationRate] = maskedClusterQualityMatlab(spike_clusters,spike_templates,pc_features,pc_feature_ind)

% based off main loop in maskedClusterQualityKilosort

% in my implementation, calling data from Kilosort/Phy is done before
% maskedClusterQuality is used

fprintf('building features matrix from clusters/templates\n')
spike_clusters = spike_clusters + 1; % because in Python indexes start at 0

% now we have to construct a new pc_features that has the features
% arranged according to *cluster* rather than template.
spike_templates = spike_templates + 1;

% pc_feature_ind has channel indexes of each template and is zero-indexed due to python
pc_feature_ind = pc_feature_ind + 1;


clusterIDs = unique(spike_clusters);
nClusters = length(clusterIDs);
nSpikes = length(spike_clusters);
nFet = 4; nFetPerChan = size(pc_features,2);
nTemplates = size(pc_feature_ind,1);

% use nFet features to calculate mahalanobis distance
newFet = zeros(nSpikes, nFetPerChan, nFet);
newFetInds = zeros(nClusters, nFet);
%tempNums = 1:nTemplates;

for c = 1:length(clusterIDs)
    %         fprintf(1, '%d/%d\n', c, length(clusterIDs))
    thisID = clusterIDs(c);

    theseSpikes = spike_clusters==thisID;
    theseTemplates = spike_templates(theseSpikes); % how many spike templates are included in these clusters?
    [inclTemps, inst] = countUnique(theseTemplates); % count how many times each template is used in cluster

    thisTemplate = inclTemps(inst==max(inst),1); % use most common template in cluster

    theseChans = pc_feature_ind(thisTemplate,1:nFet); % look at the channels with highest responses for this template

    newFetInds(c,:) = theseChans; % newFetInds stores which channels to look at for each cluster

    %subPCFetInd = pc_features(theseSpikes,:,:);


    for f = 1:nFet
        % find PCs that use the same channel
        thisChanInds = pc_feature_ind==theseChans(f);
        [chanInds,tempsWithThisChan] = find(thisChanInds');
        %spikesWithThisFet = ismember(theseTemplates, tempsWithThisChan);

        % use only templates that 1) are part of the same cluster, and 2)
        % share the same channel as thisTemplate (most common one in
        % cluster)
        inclTempsWithThisFet = find(ismember(inclTemps, tempsWithThisChan));

        % for each template that shares a channel with any of the templates
        % within the cluster of interest, including itself(?)
        for t = 1:numel(inclTempsWithThisFet)
            thisSubTemp = inclTemps(inclTempsWithThisFet(t));
            thisTfetInd = chanInds(tempsWithThisChan==thisSubTemp);

            % pc_features size is [nSpks x nFeatPerChans x nFeats], with the nFeatPerChans top
            % features projected onto channels specified by
            % pc_features_ind [nTemplates,nFeats]
            newFet(theseSpikes&spike_templates==thisSubTemp,:,f) = ...
                pc_features(theseSpikes&spike_templates==thisSubTemp,:,thisTfetInd);
        end


    end
end

pc_features = newFet;
pc_feature_ind = newFetInds;

assert(numel(size(pc_features)) == 3)

fprintf(1, 'computing cluster qualities...\n');
[clusterIDs, unitQuality, contaminationRate] = maskedClusterQualitySparse(spike_clusters, pc_features, pc_feature_ind);



