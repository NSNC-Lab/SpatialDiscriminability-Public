function [ctrl_data,laser_data,ctrl_folder,laser_folder] = ...
    sortSpikesbyTrial(database,subject,TMR,power,targetdB,spkfile)

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

subfolder = fullfile(database,subject);

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

if isfile(ctrl_data), return; end

% Load timestamps for trials

load(fullfile(subfolder,strcat(subject,'_',num2str(power),'mW_',num2str(targetdB),'dBtarget_',num2str(TMR),'dBTMR_ts.mat')));

% Load spike data
in = load(spkfile,'Mua','snips','codes','labels');

% Remove remaining laser-artifact spikes from data
[Mua,snips,codes,labels] = denoiseKsortSpks(in.Mua,in.snips,in.codes,in.labels);

clearvars in

Chans = find(~cellfun(@isempty,Mua));

%% Organize spikes by trial

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
