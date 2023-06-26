function [clustfile] = fetchSpikeClustMetrics(phy_folder)

addpath(genpath('spikes-master'));
addpath(genpath('npy-matlab'));

ksDir = fullfile(phy_folder);

params = struct;
% params.loadPCs = false;
params.loadPCs = true;

% cluster_info.tsv has info on which channel has peak amplitude per cluster
% and manual labels from Phy (good / mua / noise)

fid = fopen(fullfile(ksDir,'cluster_info.tsv'));
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s','delimiter','\t','emptyvalue',0);
fclose(fid);

[cids2, uQ, cR] = sqKilosort.maskedClusterQuality(ksDir);

ISIviol = sqKilosort.isiViolations(ksDir);

% one last check to see if units are actually MUA or single-unit (based on
% % of ISIs below a certain value)

% fetch information from good and MUA clusters
inds = find(~strcmp(C{9},'noise') & ~cellfun(@isempty,C{9})); inds(1) = [];

cids = str2double(C{1}(inds));  % good cluster IDs % starts from 0 due to python's indexing

%cids2 starts from 1
IsolationDists = uQ(ismember(cids2-1,cids));
ContamRatios = cR(ismember(cids2-1,cids));
ISIviol = ISIviol(ismember(cids2-1,cids));
IDs = cids;
chs = str2double(C{6}(inds)) + 1; % channel # for each cluster ID (need to add the 1 because Phy counts from channel 0)

clustfile = fullfile(phy_folder,'clusterMetrics.mat');

save(clustfile,'IsolationDists','ContamRatios','ISIviol','IDs','chs');

end