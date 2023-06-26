function wf = getWaveForms_JN(gwfparams,bpsignal)
% function wf = getWaveForms_JN(gwfparams,bpsignal)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
%
%%%%%%%%%
%
% Jio's edits:
% - removed spike times where the endpoint of the waveform would be
%   past the end of the recording
% - no longer calculates the mean waveform for each cluster at each channel
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));

chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

% Read spike time-centered waveforms
unitIDs = unique(gwfparams.spikeClusters);
numUnits = size(unitIDs,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);

for curUnitInd = 1:numUnits
    curUnitID = unitIDs(curUnitInd);
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
        
    % remove spikeTimeKeeps + wfWin(end) > size(bpsignal,2) - JN
    curSpikeTimes(curSpikeTimes + gwfparams.wfWin(end) > size(bpsignal,2)) = [];

    curUnitnSpikes = size(curSpikeTimes,1);
    spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
    spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
    for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
        tmpWf = bpsignal(1:gwfparams.nCh,spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(1):spikeTimeKeeps(curUnitInd,curSpikeTime)+gwfparams.wfWin(end));
        waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
    end
    % waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
    % disp(['Completed ' int2str(curUnitInd) ' units of ' int2str(numUnits) '.']);
end

% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
% wf.waveFormsMean = waveFormsMean;

end
