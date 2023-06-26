
addpath('cSPIKE');
addpath(genpath('subfunctions'));
InitializecSPIKE;
addpath(genpath('npy-matlab'));
addpath(genpath('sortingQuality-master'));

if ~isfile(fullfile(ksort_folder,'Kilosort_spks.mat'))
    spkfile = convertPhy2mat(ksort_folder);
else
    spkfile = fullfile(ksort_folder,'Kilosort_spks.mat');
end

%% Extract spikes from Kilosort pipeline and sort by trial type

% ctrl_data and laser_data are the filenames for the control and laser
% trials

% ctrl_folder and laser_folder are the folders where the matfiles and
% processed results for each condition will be stored

[ctrl_data,laser_data,ctrl_folder,laser_folder] = ...
    sortSpikesbyTrial(analysis_storage,subject,TMR,power,targetdB,spkfile);

%% Calculate performance metrics

% call processed spike files

ctrl_spks = load(ctrl_data);
laser_spks = load(laser_data);

% control
[~,goodchans,ctrl_perf] = calc_performance_SPIKE(ctrl_spks.Spks_clean,ctrl_spks.Spks_masked,...
    ctrl_folder,subject,'',[],TMR,power);

calc_performance(ctrl_spks.Spks_clean,ctrl_spks.Spks_masked,...
    ctrl_folder,subject,'',[],TMR,power); % VR-based performance

calc_performance_ISI(ctrl_spks.Spks_clean,ctrl_spks.Spks_masked,...
    ctrl_folder,subject,'',[],TMR,power);

calc_performance_RISPIKE(ctrl_spks.Spks_clean,ctrl_spks.Spks_masked,...
    ctrl_folder,subject,'',[],TMR,power);

calc_performance_spkct(ctrl_spks.Spks_clean,ctrl_spks.Spks_masked,...
    ctrl_folder,subject,'',[],TMR,power);

% laser
[~,~,laser_perf] = calc_performance_SPIKE(laser_spks.Spks_clean,laser_spks.Spks_masked,...
    laser_folder,subject,'laser',[],TMR,power);

calc_performance(laser_spks.Spks_clean,laser_spks.Spks_masked,...
    laser_folder,subject,'laser',[],TMR,power); % VR-based performance

calc_performance_ISI(laser_spks.Spks_clean,laser_spks.Spks_masked,...
    laser_folder,subject,'laser',[],TMR,power);

calc_performance_RISPIKE(laser_spks.Spks_clean,laser_spks.Spks_masked,...
    laser_folder,subject,'laser',[],TMR,power);

calc_performance_spkct(laser_spks.Spks_clean,laser_spks.Spks_masked,...
    laser_folder,subject,'laser',[],TMR,power);

%% Plot spatial grids of all channels

subfolder = fullfile(analysis_storage,subject);

close all;
plotGrids(ctrl_perf,laser_perf,TMR,power,subfolder)

%% Plot rasters from hotspots only

Chans = findResponsiveChans(ctrl_spks);

response_folder = [subfolder filesep 'Hotspot responses'];
if ~isfolder(response_folder)
    mkdir(response_folder); 
end

for ch = goodchans
    [inds_clean] = findHotspots(ctrl_perf.Max{ch},ctrl_perf.d{ch},0,70,1);
    [inds_masked] = findHotspots(ctrl_perf.max_masked{ch},ctrl_perf.d_masked{ch},0,70,1);

    if ~isempty(inds_clean)
        for ii = inds_clean'
        perfs = [ctrl_perf.Max{ch}(ii) laser_perf.Max{ch}(ii)];
        lind = ii;

        [title_str] = plotCompareResponses(ctrl_spks.Spks_clean{ch}(:,ii,:),laser_spks.Spks_clean{ch}(:,ii,:),0.5,3.5,ch,lind,perfs);
        saveas(gcf,[response_folder filesep title_str '.png']);
        close;
        end
    end

    if ~isempty(inds_masked)
        for ii = inds_masked'

        [mloc,tloc] = ind2sub([4 4],ii);
        mloc = 5 - mloc;

        perfs = [ctrl_perf.max_masked{ch}(ii) laser_perf.max_masked{ch}(ii)];
        lind = [tloc mloc];

        [title_str] = plotCompareResponses(ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:),...
            laser_spks.Spks_masked{ch}(:,tloc,mloc,:),0.5,3.5,ch,lind,perfs);
        saveas(gcf,[response_folder filesep title_str '.png']);
        close;
        end
    end
end

%% Determine single and multi-units using sortingQuality toolbox

[cids,vISI,uQ,cR] = runSpikeClusterMetrics(ksort_folder,Chans);

% for clusters where mahalanobis distance can't be calculated, use ISI
% violations as only criteria for single unit activity
uQ(uQ == 0) = 50; cR(isnan(cR)) = 0.01;

% % relaxed thresholds
SU = uQ > 15 & cR < 0.25 & vISI' < 0.05;

% Narrow spiking or regular spiking?
spon_FR = []; spk_width = []; 

for ch = Chans

    spk_width = cat(1,spk_width,calcSpkWidth(ctrl_spks.snips_clean{ch},ctrl_spks.snips_masked{ch}));

    temp = [];
    for t = 1:4
        temp = cat(1,temp,squeeze(calcFRPeriod(ctrl_spks.Spks_clean{ch}(:,t,:),'spontaneous')));
    end
    for t = 1:16
        [tloc,mloc] = ind2sub([4 4],t);
        temp = cat(1,temp,squeeze(calcFRPeriod(ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:),'spontaneous')));
    end
    spon_FR = cat(1,spon_FR,mean(temp));
end
    
ns_inds = spk_width < 0.5;

% table
unitInfo = array2table([cids,uQ,cR,vISI',SU,ns_inds],'VariableNames',{'CH','IsolDist','ContamRatio','ISIViol','SingleUnit?','NarrowSpk?'})
save([subfolder filesep 'unitInfo.mat'],'unitInfo');
