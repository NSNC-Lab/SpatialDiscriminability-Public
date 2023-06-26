function [spkfile] = convertPhy2mat(phy_folder)

% use code to convert phy data into spike times for discriminability
% analysis

addpath(genpath(fullfile('..','spikes-master')));
addpath(genpath(fullfile('..','npy-matlab')));

ksDir = fullfile(phy_folder);

params = struct;
% params.loadPCs = false;
params.loadPCs = true;

% spikeStruct includes information for all non-noise spikes
spikeStruct = loadKSdir(ksDir,params);

% cluster_info.tsv has info on which channel has peak amplitude per cluster
% and manual labels from Phy (good / mua / noise)

fid = fopen(fullfile(ksDir,'cluster_info.tsv'));
C = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s','delimiter','\t','emptyvalue',0);
fclose(fid);

winlen = 64;

[cids2, uQ, cR] = sqKilosort.maskedClusterQuality(phy_folder);

gwfparams.fileName = 'Raws.i16';

gwfparams.dataDir = phy_folder;    % KiloSort/Phy output folder
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-round(winlen/4)+1 winlen-round(winlen/4)];                % Number of samples before and after spiketime to include in waveform

fileName = fullfile(gwfparams.dataDir,gwfparams.fileName);           
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel

% one last check to see if units are actually MUA or single-unit (based on
% % of ISIs below a certain value)

% fetch information from good and MUA clusters
inds = find(~strcmp(C{9},'noise') & ~cellfun(@isempty,C{9})); inds(1) = [];

cids = str2double(C{1}(inds));  % good cluster IDs % starts from 0 due to python's indexing

%cids2 starts from 1
IsolationDists = uQ(ismember(cids2-1,cids));
ContamRatios = cR(ismember(cids2-1,cids));

% for clu = 1:length(cids)
%     C{9}{inds(clu)} = checkifMua(spikeStruct.st(spikeStruct.clu == cids(clu)),labels{clu});
% end

chs = str2double(C{6}(inds)) + 1; % channel # for each cluster ID (need to add the 1 because Phy counts from channel 0)
n_spikes = str2double(C{10}(inds));

% load band-passed signal outside of getWaveforms to save space on
% waveForms matrix

disp('Loading and filtering raw signal...');

tic;

fID = fopen(fileName);
x = fread(fID,[gwfparams.nCh nSamp],'int16');
fclose(fID);

% create butterworth band-pass filter based on synapse parameters
[b,a] = butter(3, [300 5000] / (spikeStruct.sample_rate/2) );
for ch = 1:gwfparams.nCh
    x(ch,:) = filtfilt(b,a,x(ch,:));
end

elapsedTime = toc;
disp(['Elapsed raw-signal loading and filtering time was ' num2str(elapsedTime) ' s.']);

% arrange spike times to channel by cluster ID
Mua = cell(1,32);
snips = cell(1,32);
codes = cell(1,32); % cluster ID (from 0)
labels = cell(1,32); % from SU or MU?

chMap = double(readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy')) + 1);

for ci = 1:length(cids)
            
    gwfparams.nWf = n_spikes(ci);                    % Number of waveforms per unit to pull out
    gwfparams.spikeTimes = round(spikeStruct.st(ismember(spikeStruct.clu,cids(ci)))*spikeStruct.sample_rate); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = spikeStruct.clu(ismember(spikeStruct.clu,cids(ci))); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
    
    % extract waveforms for all non-noise clusters
    wf = getWaveForms_JN(gwfparams,x);
    
    % getWaveforms shuffles channels by the channel map, which is based on how
    % the data is streamed to the disk. Since we converted the TDT SEV files
    % into binary ourselves, the data is streamed as [1 2 3 ... 32]. However,
    % Kilosort discards bad channels so length(chMap) < 32
    
    % convert waveforms from int16 to double and scale snippets down by
    % 1E06, since binary file was made using that scale factor (see
    % convertRaw2Binary)
    
    % get all snippets from cluster and channel, then
    temp_snips = squeeze(wf.waveForms(:,:,chMap == chs(ci),:)) / 1e6;
    
    temp = wf.spikeTimeKeeps' / spikeStruct.sample_rate;
    
    Mua{1,chs(ci)} = cat(1,Mua{1,chs(ci)},temp);
    snips{1,chs(ci)} = cat(1,snips{1,chs(ci)},temp_snips);
    codes{1,chs(ci)} = cat(1,codes{1,chs(ci)},double(wf.unitIDs)*ones(size(temp)));
    labels{1,chs(ci)} = cat(1,labels{1,chs(ci)},repmat(C{9}(inds(ci)),size(temp_snips,1),1));
    disp(['Completed ' int2str(ci) ' units of ' int2str(length(cids)) '.']);
end

% sort spiketimes in time for non-empty channels

goodchans = find(~cellfun(@isempty,Mua));

for ch = 1:length(goodchans)
    [Mua{1,goodchans(ch)},time_inds] = sort(Mua{1,goodchans(ch)});
    snips{1,goodchans(ch)} = snips{goodchans(ch)}(time_inds,:);
    codes{1,goodchans(ch)} = codes{goodchans(ch)}(time_inds);
    labels{1,goodchans(ch)} = labels{goodchans(ch)}(time_inds);
end

spkfile = fullfile(phy_folder,'Kilosort_spks.mat');

save(spkfile,'snips','Mua','codes','labels');

% next: run the script that sorts the spikes and snippets by trial

end