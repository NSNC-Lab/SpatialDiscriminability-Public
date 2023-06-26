function [ctrl_data,laser_data,ctrl_folder,laser_folder] = ...
    extract_spatial_grid_Kilosort_V3(database,tanks,subject,TMR,power,targetdB,spkfile)

% ctrl_data and laser_data are matfiles that contain:

% Spks: spike times for each trial
% n_cum: PSTH for spatial grid location
% num_spks: # spikes in each trial per location (spontaneous + driven)
% snips: spike waveforms for each spike in Spk
% clu_IDs: cluster ID from Kilosort
% labels: spk comes from good (single-unit) or MUA, based on 1) Kilosort
% and 2) whether enough ISIs in the cluster are below a certain value

t_before = 1;
t_after = 4;
t_bin = 1/50;
t_vec = -t_before:t_bin:t_after;

subfolder = fullfile(database,subject,[subject,'-',num2str(power),'mW_',num2str(targetdB),'dBtarget-',num2str(TMR),'dBTMR_Kilosort_correctMap']);
if ~isfolder(subfolder)
    mkdir(subfolder)
end

ctrl_folder = [subfolder filesep 'No Laser'];
if ~isfolder(ctrl_folder)
    mkdir(ctrl_folder)
end

laser_folder = [subfolder filesep 'Laser'];
if ~isfolder(laser_folder)
    mkdir(laser_folder)
end

ctrl_data = strcat(ctrl_folder,'/',subject,'_',num2str(power),'mW_',num2str(targetdB),'dBtarget_',num2str(TMR),'dBTMR_raw(-1,4)_Kilosort.mat');
laser_data = strcat(laser_folder,'/',subject,'_',num2str(power),'mW_',num2str(targetdB),'dBtarget_',num2str(TMR),'dBTMR_rawlaser(-1,4)_Kilosort.mat');

if isfile(ctrl_data) && isfile(laser_data), return; end

clearvars TSClean TSMasked TSMaskedlaser TSCleanlaser

TSClean{4,2}=[];
TSMasked{4,4,2}=[];
TSCleanlaser{4,2}=[];
TSMaskedlaser{4,4,2}=[];

% since the spkfile treats the blocks as if the experiment were a
% continuous whole, we have to account for the length of time each block
% adds to the experiment/timestamps
exp_time = 0;

for nT = 1:numel(tanks)
    
    % import TM info
    data = TDTbin2mat(tanks{nT},'TYPE',5);   % extract scalars only
    mloc = data.scalars.MskS.data + 1;
    mt = data.scalars.MskV.data;
    tloc = data.scalars.TarS.data + 1;
    tt = data.scalars.TarV.data;
    
    % NOTE: TARGET/MASKER LOCATIONS CHANGE ONE SECOND BEFORE STIMULUS PLAYS, SO
    % ADD T_BEFORE TO ALL TS VALUES
        
    TS = data.scalars.TarV.ts + t_before + exp_time;
    
    data_stream = TDTbin2mat(tanks{nT},'TYPE',4,'STORE','Raws','Channel',1); 

    exp_time = exp_time + size(data_stream.streams.Raws.data,2)/24414.0625;  % update exp_time
    
    % Sometimes the TSes overcount by one
    trial_lengths = round(diff(TS));
    
    % For some reason, the TSes for stimulus presentation will have a dummy
    % trial where mt == 31 or mt == 0, which is impossible given the logic of the circuit.
    % This doesn't affect the experimental data thankfully, but we need to take
    % it out of the TSes anyway.

    % if a trial is erroneously added (usually happens right at the first
    % trial(=)
    if sum(trial_lengths - 5) > 0
        bad = find(trial_lengths ~= 5,1) + 1;
        mloc(bad) = [];
        tloc(bad) = [];
        mt(bad) = [];
        tt(bad) = [];
        TS(bad) = [];
    end
    
    % if an indexing thing goes wrong mid-way through the recording
    badind = find(tt == 31 | tt == 0);
    if ~isempty(badind)
        mloc(badind) = [];
        tloc(badind) = [];
        mt(badind) = [];
        tt(badind) = [];
        TS(badind) = [];
    end
    
    %% Find timestamps for each stimuli/TM location trial
    
    % for partial blocks or blocks with extra trials
    nTS = length(TS);
    
    % Use location mat to organize timestamps by TM and laser
    laser = TDTbin2mat(tanks{nT},'TYPE',2,'STORE','Buff');

    % the first index is a dummy entry for cueing the correct laser
    % trial order
    laserinfo = laser.epocs.Buff.data(2:end);
    
    % if there is a bad trial index after removing erroneously added
    % trials, remove it from laser info so that indexing isn't messed up
    if length(laserinfo) > nTS && ~isempty(badind)
        laserinfo(badind) = [];
    elseif length(laserinfo) > nTS %% if these lengths are not equal AND there's no badind, I have no idea what's up
        minlen = min([length(laserinfo) nTS+1]);
        laserinfo = laserinfo(1:minlen);
    end
    
    if ( (rem(length(TS),80) ~= 0 || rem(length(laserinfo),80) ~= 0)) && (length(TS) > 80 && length(TS) < 320) || ...
            length(TS) > 400 
        nTS = 80*floor(length(TS)/80);
        
        mloc(nTS+1:end) = [];
        tloc(nTS+1:end) = [];
        mt(nTS+1:end) = [];
        tt(nTS+1:end) = [];
        TS(nTS+1:end) = [];
        laserinfo(nTS+1:end) = [];
    end
        
    % if the block almost has a full set of trials, that's fine, but take out
    % the last time stamp which might be a partial trial
    
    % check for unique TS
    alltrials = [tloc',mloc',tt',mt',laserinfo(1:nTS)];
    uniquetrials = unique(alltrials,'rows');
    if rem(size(uniquetrials,1) - 16,64) ~= 0
        warning('Amount of trials per type is not equal.');
    end
    
    % each block should have a total of 5 trials per trial type. sometimes,
    % synapse prematurely stops recording the experiment. for blocks that are
    % just close to either having 1 or 5 trial per trial type, we'll discard the
    % last trial
    
%     % if the amount of trials in the block is almost 400 or almost 80
%     if (length(TS) >= 390 || length(TS) < 80) && length(TS) < 400
%         nTS = nTS - 1;
%     end
    
    for n = 1:nTS
        if laserinfo(n) == 0
            for j = 1:2
                if (tt(n) == j) && (mt(n) == 0) && (mloc(n) == 16)   % clean trial
                    TSClean{tloc(n),j} = cat(1,TSClean{tloc(n),j},TS(n));
                elseif (tt(n) == j) && (mt(n) ~= 0)   % masked trials, diff locations
                    TSMasked{tloc(n),mloc(n),j} = cat(1,TSMasked{tloc(n),mloc(n),j},TS(n));
                elseif tt(n) >= 10*j+1 && tt(n) <= 10*(j+1)     % colocated
                    TSMasked{tloc(n),tloc(n),j} = cat(1,TSMasked{tloc(n),tloc(n),j},TS(n));
                end
            end
        elseif laserinfo(n) == 1
            for j = 1:2
                if (tt(n) == j) && (mt(n) == 0) && (mloc(n) == 16)
                    TSCleanlaser{tloc(n),j} = cat(1,TSCleanlaser{tloc(n),j},TS(n));
                elseif (tt(n) == j) && (mt(n) ~= 0)
                    TSMaskedlaser{tloc(n),mloc(n),j} = cat(1,TSMaskedlaser{tloc(n),mloc(n),j},TS(n));
                elseif tt(n) >= 10*j+1 && tt(n) <= 10*(j+1)
                    TSMaskedlaser{tloc(n),tloc(n),j} = cat(1,TSMaskedlaser{tloc(n),tloc(n),j},TS(n));
                end
            end
        end
    end
    
    trialinfo{nT} = alltrials;
    
    clearvars data

end

ts_data = fullfile(subfolder,...
    strcat(subject,'_',num2str(power),'mW_',num2str(targetdB),'dBtarget_',num2str(TMR),'dBTMR_ts'));
save(ts_data,'trialinfo');

% skip function if matfiles are already made

% if isfile(ctrl_data) && isfile(laser_data), return; end

% Load spike data
in = load(spkfile,'Mua','snips','codes','labels');

% Remove remaining laser-artifact spikes from data
[Mua,snips,codes,labels] = denoiseKsortSpks(in.Mua,in.snips,in.codes,in.labels);

clearvars in

Chans = find(~cellfun(@isempty,Mua));
%% Organize LFPs and spikes by trial

n_cum_clean = cell(1,32);
n_cum_masked = cell(1,32);
n_cum_cleanlaser = cell(1,32);
n_cum_maskedlaser = cell(1,32);

num_spks_clean = cell(1,32);
num_spks_masked = cell(1,32);
num_spks_cleanlaser = cell(1,32);
num_spks_maskedlaser = cell(1,32);

Spks_clean = cell(1,32);
Spks_masked = cell(1,32);
Spks_cleanlaser = cell(1,32);
Spks_maskedlaser = cell(1,32);

snips_clean = cell(1,32);
snips_masked = cell(1,32);
snips_cleanlaser = cell(1,32);
snips_maskedlaser = cell(1,32);

clu_IDs_clean = cell(1,32);
clu_IDs_masked = cell(1,32);
clu_IDs_cleanlaser = cell(1,32);
clu_IDs_maskedlaser = cell(1,32);

for ch = Chans
    
    % clean data
    
    for x = 1:4
        for z = 1:2
            
            [raster_clean{1,ch}{x,z},n_cum_clean{1,ch}{x,z}]=...
                cal_raster(TSClean{x,z},Mua{1,ch},t_vec);
            [raster_cleanlaser{1,ch}{x,z},n_cum_cleanlaser{1,ch}{x,z}]=...
                cal_raster(TSCleanlaser{x,z},Mua{1,ch},t_vec);
            
            [waveforms_clean{1,ch}{x,z},codes_clean{1,ch}{x,z},spktype_clean{1,ch}{x,z}] = sort_snippets(...
                TSClean{x,z},Mua{1,ch},snips{1,ch},codes{1,ch},labels{1,ch},t_vec);
            [waveforms_cleanlaser{1,ch}{x,z},codes_cleanlaser{1,ch}{x,z},spktype_cleanlaser{1,ch}{x,z}] = sort_snippets(...
                TSCleanlaser{x,z},Mua{1,ch},snips{1,ch},codes{1,ch},labels{1,ch},t_vec);
            
            for k=1:length(raster_clean{1,ch}{x,z})
                Spks_clean{1,ch}{k,x,z} = raster_clean{1,ch}{x,z}{1,k};
                num_spks_clean{1,ch}{x,z}(k) = length(raster_clean{1,ch}{x,z}{1,k});
                snips_clean{1,ch}{k,x,z} = waveforms_clean{1,ch}{x,z}{1,k};
                clu_IDs_clean{1,ch}{k,x,z} = codes_clean{1,ch}{x,z}{1,k};
                labels_clean{1,ch}{k,x,z} = spktype_clean{1,ch}{x,z}{1,k};
            end
            
            for k=1:length(raster_cleanlaser{1,ch}{x,z})
                Spks_cleanlaser{1,ch}{k,x,z}=raster_cleanlaser{1,ch}{x,z}{1,k};
                num_spks_cleanlaser{1,ch}{x,z}(k)=length(raster_cleanlaser{1,ch}{x,z}{1,k});
                snips_cleanlaser{1,ch}{k,x,z} = waveforms_cleanlaser{1,ch}{x,z}{1,k};
                clu_IDs_cleanlaser{1,ch}{k,x,z} = codes_cleanlaser{1,ch}{x,z}{1,k};
                labels_cleanlaser{1,ch}{k,x,z} = spktype_cleanlaser{1,ch}{x,z}{1,k};
            end
        end
    end
    
    % mixed data
    
    if ~isempty(TSMasked{1,1,1})
        
        for x = 1:4 % target location (90 -> -90)
            for y = 1:4 % masker location  (90 -> -90)
                for z = 1:2 % target identity 
                    
                    [raster_masked{1,ch}{x,y,z},n_cum_masked{1,ch}{x,y,z}]=...
                        cal_raster(TSMasked{x,y,z},Mua{1,ch},t_vec);
                    [raster_maskedlaser{1,ch}{x,y,z},n_cum_maskedlaser{1,ch}{x,y,z}]=...
                        cal_raster(TSMaskedlaser{x,y,z},Mua{1,ch},t_vec);
                    
                    [waveforms_masked{1,ch}{x,y,z},codes_masked{1,ch}{x,y,z},spktype_masked{1,ch}{x,y,z}] = sort_snippets(...
                        TSMasked{x,y,z},Mua{1,ch},snips{1,ch},codes{1,ch},labels{1,ch},t_vec);
                    [waveforms_maskedlaser{1,ch}{x,y,z},codes_maskedlaser{1,ch}{x,y,z},spktype_maskedlaser{1,ch}{x,y,z}] = sort_snippets(...
                        TSMaskedlaser{x,y,z},Mua{1,ch},snips{1,ch},codes{1,ch},labels{1,ch},t_vec);
                    
                    for k=1:length(raster_masked{1,ch}{x,y,z})
                        Spks_masked{1,ch}{k,x,y,z}=raster_masked{1,ch}{x,y,z}{1,k};
                        num_spks_masked{1,ch}{x,y,z}(k)=length(raster_masked{1,ch}{x,y,z}{1,k});
                        snips_masked{1,ch}{k,x,y,z} = waveforms_masked{1,ch}{x,y,z}{1,k};
                        clu_IDs_masked{1,ch}{k,x,y,z} = codes_masked{1,ch}{x,y,z}{1,k};
                        labels_masked{1,ch}{k,x,y,z} = spktype_masked{1,ch}{x,y,z}{1,k};
                    end

                    for k=1:length(raster_maskedlaser{1,ch}{x,y,z})
                        Spks_maskedlaser{1,ch}{k,x,y,z}=raster_maskedlaser{1,ch}{x,y,z}{1,k};
                        num_spks_maskedlaser{1,ch}{x,y,z}(k)=length(raster_maskedlaser{1,ch}{x,y,z}{1,k});
                        snips_maskedlaser{1,ch}{k,x,y,z} = waveforms_maskedlaser{1,ch}{x,y,z}{1,k};
                        clu_IDs_maskedlaser{1,ch}{k,x,y,z} = codes_maskedlaser{1,ch}{x,y,z}{1,k};
                        labels_maskedlaser{1,ch}{k,x,y,z} = spktype_maskedlaser{1,ch}{x,y,z}{1,k};
                    end

                end
            end
        end
        
    end
    
end


save(ctrl_data,'n_cum_clean','Spks_clean','num_spks_clean','n_cum_masked','Spks_masked','num_spks_masked',...
    'snips_clean','snips_masked','clu_IDs_clean','clu_IDs_masked','labels_clean','labels_masked');

n_cum_clean = n_cum_cleanlaser;
Spks_clean = Spks_cleanlaser;
n_cum_masked = n_cum_maskedlaser;
Spks_masked = Spks_maskedlaser;
num_spks_clean = num_spks_cleanlaser;
num_spks_masked = num_spks_maskedlaser;

snips_clean = snips_cleanlaser;
snips_masked = snips_maskedlaser;
clu_IDs_clean = clu_IDs_cleanlaser;
clu_IDs_masked = clu_IDs_maskedlaser;
labels_clean = labels_cleanlaser;
labels_masked = labels_maskedlaser;

save(laser_data,'n_cum_clean','Spks_clean','num_spks_clean','n_cum_masked','Spks_masked','num_spks_masked',...
    'snips_clean','snips_masked','clu_IDs_clean','clu_IDs_masked','labels_clean','labels_masked');

end
