function [cids,ISI_viol,uQ,cR] = runSpikeClusterMetrics(ksort_folder,good_chs)

% input: ksDir (kilosort data stored here)

addpath(genpath('spikes-master'));
addpath(genpath('npy-matlab'));

params = struct;
params.loadPCs = true;

% spikeStruct includes information for all non-noise spikes
spikeStruct = loadKSdir(ksort_folder,params);
spkfile = fullfile(ksort_folder,'Kilosort_spks.mat');
% Load spike data
in = load(spkfile,'Mua','snips','codes','labels');

% Remove remaining laser-artifact spikes from data
[Mua,~,codes,~] = denoiseKsortSpks(in.Mua,in.snips,in.codes,in.labels);

spikeTimes = round(spikeStruct.st*24414.0625)/24414.0625; % Vector of cluster spike times (in samples) same length as .spikeClusters
spikeClusters = double(spikeStruct.clu); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

bad_chans = setdiff(1:32,good_chs);
for b = bad_chans
    Mua{b} = [];
    codes{b} = [];
end

% sort all remaining spikes from Mua
[Mua_new,inds] = sort(vertcat(Mua{:}));
temp = vertcat(codes{:}); codes_temp = temp(inds);

% only include clusters in both spikeStruct and codes

% sometimes the saved spikes have a cluster not included in kilosort
% anymore/labeled as noise, need to investigate this

real_clusters = intersect(codes_temp,spikeClusters);
ind_codes = ismember(codes_temp,real_clusters); % part of clusters we wanna keep

Mua_new = Mua_new(ind_codes); codes_new = codes_temp(ind_codes);

chs = [];
for c = 1:32
    chs = cat(1,chs,c*ones(numel(Mua{c}),1));
end
chs = chs(inds); chs(~ind_codes) = []; 

% find spikes from Mua_new in spikeTimes (old, not denoised)
inds_grp = ismember([spikeTimes,spikeClusters],[Mua_new,codes_new],'rows');

real_spks = [spikeTimes(inds_grp),spikeClusters(inds_grp)];
PCs = spikeStruct.pcFeat(inds_grp,:,:);
spikeTemplates = spikeStruct.spikeTemplates(inds_grp);

% after grouping clusters in same channel
[cids,uQ,cR] = maskedClusterQualityMatlab(chs-1,spikeTemplates,PCs,spikeStruct.pcFeatInd);
ISI_viol = [];

nonemps = find(~cellfun(@isempty,Mua));
for s = 1:length(nonemps)
    ISIs = diff(Mua{nonemps(s)});
    ISI_viol(s) = sum(ISIs < 0.002)/numel(ISIs);
end

end