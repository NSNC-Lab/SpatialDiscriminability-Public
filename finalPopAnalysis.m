
% must have run:
% 1) spgrid_automate_KSort
% 2) spgrid_loop_analysis
% 3) spgrid_LFP_analysis for each subject

clear; close all; clc

addpath('stimuli');
addpath(genpath('subfunctions'));
addpath('sigstar-master');
addpath(genpath('sortingQuality-master'));

analysis_storage = 'Arch-data/analysis'; %uigetdir(cd,'Select where to save processed spikes and results.');   % 
data_storage = 'Arch-data/processed-Kilosort'; %uigetdir(cd,'Select folder with stored data tanks.');     

addpath(analysis_storage);

subjects = {'616283','616299','616302','616303','616304','616305','616306','616307','616312'};

perf_thresh = 70;
d_thresh = 1;
whichpower = 10;
whichTMR = 0;

exps = struct;
perf_methods = {'','_SPIKE','_RI_SPIKE','_ISI','_spkct'};

for ns = 1:length(subjects)

    subfolder = dir(fullfile(analysis_storage,subjects{ns}));

    expind = find(contains({subfolder.name},[num2str(whichpower) 'mW']) & ...
        contains({subfolder.name},[num2str(whichTMR) 'dBTMR']) & ...
        contains({subfolder.name},'Kilosort_correctMap'));

    ctrl_folder = fullfile(subfolder(expind).folder,subfolder(expind).name,'No Laser');
    laser_folder = fullfile(subfolder(expind).folder,subfolder(expind).name,'Laser');

    addpath(ctrl_folder);
    addpath(laser_folder);

    ctrl_files = dir(ctrl_folder);
    spks_file = ctrl_files(contains({ctrl_files.name},'(-1,4)')).name;
    exps(ns).ctrl_spks = load([ctrl_folder filesep spks_file],'Spks_clean','Spks_masked','num_spks_clean','num_spks_masked','n_cum_clean','n_cum_masked','snips_clean','snips_masked');

    for pf = 1:length(perf_methods)
        perf_file = ctrl_files(contains({ctrl_files.name},['TMR_performance' perf_methods{pf} '.mat'])).name;
        exps(ns).(['ctrl_perf' perf_methods{pf}]) = load([ctrl_folder filesep perf_file],'Max','max_masked','d','d_masked','clean_performance','masked_performance','null_clean','null_masked','opt_tau_masked','opt_tau');
    end

    laser_files = dir(laser_folder);
    spks_file = laser_files(contains({laser_files.name},'(-1,4)')).name;
    exps(ns).laser_spks = load([laser_folder filesep spks_file],'Spks_clean','Spks_masked','num_spks_clean','num_spks_masked','n_cum_clean','n_cum_masked','snips_clean','snips_masked');

    for pf = 1:length(perf_methods)
        perf_file = laser_files(contains({laser_files.name},['TMR_performance' perf_methods{pf} 'laser.mat'])).name;
        exps(ns).(['laser_perf' perf_methods{pf}]) = load([laser_folder filesep perf_file],'Max','max_masked','d','d_masked','clean_performance','masked_performance','null_clean','null_masked','opt_tau_masked','opt_tau');
    end

    trialfiles = dir(fullfile(subfolder(expind).folder,subfolder(expind).name));
    trial_mat = trialfiles(contains({trialfiles.name},'.mat')).name;
    load(fullfile(subfolder(expind).folder,subfolder(expind).name,trial_mat));
    exps(ns).trialinfo = restoreTS(vertcat(trialinfo{:}));

    csdind = find(contains({subfolder.name},'10mW') & ...
        contains({subfolder.name},'0dBTMR') & ...
        contains({subfolder.name},'LFPs'));

    csd_folder = fullfile(subfolder(csdind).folder,subfolder(csdind).name,'No Laser');
    load([csd_folder filesep 'CSD_results.mat']);

    exps(ns).gran_layers = gran_layers + 1; % add 1 because gran_layers is indexed between 1 and 6, but CSD is calculated between channels 2 and 7
    exps(ns).csd_onset = csd_onset;

    clc; clearvars trialinfo

end
fs = 24414.0625; clc

%% Find SU/MU and NS/RS channels first

artifact_chs = {[],29,22,25,[],27,[23 26],[],25};

% all cells in each response type - not used in paper
on_sus = {[2 25],[],[],[],[],[],[],[],[]};
on = {[],[],[],[],[],24,[],[],[]};
off = {[17 19],[23 26],[],[],[],[],[],[],[]};
supp = {[27],[],22,[],[],[19 20 21 28 31 32],26,24,[]};
ramp = {8,31,[],[],[],[],[],[],[]};
weak = {[22 24],[17 25 29],[20 24 28],...
    [22 25],[],[15 27],[],[],[]};
NS_chs = {[9 14 20 23 29],[3 24],[],[],[],[23 30],[21 23 24 31],[8 11],[25]};
ls_units = {18,[],27,[],[],[],[22 24 25 30],[],[]};

cellinfo = struct;
for ns = 1:length(exps)
    cellinfo(ns).narrow_spiking = NS_chs{ns};
    cellinfo(ns).reg_spiking = sort(horzcat(on_sus{ns},on{ns},off{ns},supp{ns},...
        ramp{ns},NS_chs{ns},ls_units{ns},weak{ns}));
end

layer_names = {'L1','L2/3','L4','L5','L6'};

% for kilosort metrics
if whichpower == 10 && whichTMR == 0
    for ns = 1:length(subjects)
        % look through tanks from one subject
        subject_dir = dir(fullfile(tank_storage,subjects{ns}));
        subject_dir(~contains({subject_dir.name},subjects{ns})) = [];

        % search withink tanks to see if Kilosort is saved and if the power and
        % TMR match whichpower and whichTMR
        for nT = 1:length(subject_dir)

            % look through each tank in the subject folder
            tank_dir = dir(fullfile(subject_dir(nT).folder,subject_dir(nT).name));

            power = str2double(extractBetween(subject_dir(nT).name,[subjects{ns} '-'],'mW'));
            TMR = str2double(extractBetween(subject_dir(nT).name,'dB-','dBTMR'));

            % for tank with Kilosort-merged folder, find other tanks with same
            % experiment parameters (dB TMR, laser power)
            if any(strcmp({tank_dir.name},'Kilosort-merged')) && power == whichpower && TMR == whichTMR

                % find the other tanks that have the same TMR and power. all
                % this information should be on the tank names which follow the
                % same naming convention
                tank_inds = find(contains({subject_dir.name},subject_dir(nT).name(1:end-14)));

                tank = fullfile(subject_dir(tank_inds(1)).folder,subject_dir(tank_inds(1)).name);

                ksort_folder{ns} = fullfile(tank,'Kilosort-merged');
            end
        end
    end
end

%% Calculate SU or MU based on spike cluster metrics

if ~exist('cids','var')
    for ns = 1:length(exps)
        Chans = findResponsiveChans(exps(ns).ctrl_spks);
        Chans = setdiff(Chans,artifact_chs{ns});

        [cids{ns},vISI{ns},uQ{ns},cR{ns}] = runSpikeClusterMetrics(ksort_folder{ns},Chans);

        % if mahalanobis distance can't be calculated, assume both
        % contamination ratio and isolation distances are within their
        % respective ranges for single-unit activity
        uQ{ns}(uQ{ns} == 0) = 50; 
        cR{ns}(isnan(cR{ns})) = 0.01;
    end
    save('spikeMetrics.mat','cids','vISI','uQ','cR');
end

for ns = 1:length(exps)
    SU{ns} = cids{ns}(uQ{ns} > 15 & cR{ns} < 0.25 & vISI{ns} < 0.05);
end

% 85 single units -> 3 rejected due to low trial number

folder = 'results';
if ~isfolder(folder), mkdir(folder); end

%%

[results,perf_comps] = popSearchHotspots(exps,subjects,perf_thresh,d_thresh,cellinfo,SU,artifact_chs);

% implement logic for control and emergent hotspots
d1 = find( [results.ctrl_d] >= d_thresh & [results.ctrl_perf] >= perf_thresh);  % control
d2 = find( ([results.ctrl_d] < d_thresh & [results.ctrl_perf] < perf_thresh) & ...
    ([results.laser_d] >= d_thresh & [results.laser_perf] >= perf_thresh));   % emergent

mainset = results(d1);
emergent = results(d2);

singleunits = mainset(strcmp({mainset.unit},'single'));
labels = {'multi','single'};

% results contain all hotspots, regardless of whether or not
% 1) performance was above threshold in laser but not control (d2/emergent)
% 2) unit was single or multi

% mainset and emergent should be mutually exclusive

save('mainset.mat','mainset');
save('emergent.mat','emergent');

%% S2 - RS vs NS single unit plot

plotSpikeWidthvsFR(exps,realInds);
savefig(gcf,[folder filesep 'RS and NS scatterplot.fig'])

plotNumUnitsbyLayer(exps,[subinds_mainset,vertcat(cids{:}),SU_inds2],folder);

%% Plot laser onset PSTHs for NS cells only - S5a

NS_chs = {[9 14 20],[3 24],[],[],[],23,23,[8 11],[]};

% center time vector around laser onset

dt = 2 / 1000;
% T_vec = -0.1:dt:(0+dt);
T_vec = [-fliplr(0.05:dt:0.1) (-0.05+dt):dt:(0+dt)];

% pool all NS units and determine if suppression was significant during
% first dt ms of laser

psth_ctrl = []; psth_laser = [];
for ns = 1:length(exps)   
    Chans = NS_chs{ns};
    for ch = Chans

        temp1 = cellfun(@(x) histcounts(x,T_vec),exps(ns).ctrl_spks.Spks_clean{ch},'UniformOutput',false);
        temp1 = vertcat(temp1{:});
        temp2 = cellfun(@(x) histcounts(x,T_vec),exps(ns).ctrl_spks.Spks_masked{ch},'UniformOutput',false);
        temp2 = vertcat(temp2{:});
        psth_ctrl = cat(1,psth_ctrl,mean([temp1;temp2])/dt);

        temp1 = cellfun(@(x) histcounts(x,T_vec),exps(ns).laser_spks.Spks_clean{ch},'UniformOutput',false);
        temp1 = vertcat(temp1{:});
        temp2 = cellfun(@(x) histcounts(x,T_vec),exps(ns).laser_spks.Spks_masked{ch},'UniformOutput',false);
        temp2 = vertcat(temp2{:});
        psth_laser = cat(1,psth_laser,mean([temp1;temp2])/dt);

    end
end

[ll_ctrl,ul_ctrl,av_ctrl] = sem(psth_ctrl);
[ll_laser,ul_laser,av_laser] = sem(psth_laser);

onset = T_vec == -0.05;

[h,p] = ttest2(psth_ctrl(:,onset),psth_laser(:,onset));
[d] = computeCohen_d(psth_ctrl(:,onset),psth_laser(:,onset))

figure('unit','inches','position',[5 5 2.4 1.8]);
%plot(T_vec(1:end-1)+0.05,av_ctrl/dt,'k','linewidth',1);
plot(T_vec(1:end-1)+0.05,av_ctrl,'k','linewidth',1);
%hold on; plot(T_vec(1:end-1)+0.05,av_laser/dt,'r','linewidth',1);
hold on; plot(T_vec(1:end-1)+0.05,av_laser,'r','linewidth',1);
ylabel('FR (Hz)'); xlabel('Time after laser onset (s)'); set(gca,'fontsize',8)
text(0,30,{['p = ' num2str(p)],['d = ' num2str(d)]});

patch([T_vec(1:end-1) fliplr(T_vec(1:end-1))]+0.05,[ll_ctrl fliplr(ul_ctrl)],'black','facealpha',0.3,'edgecolor','none')
patch([T_vec(1:end-1) fliplr(T_vec(1:end-1))]+0.05,[ll_laser fliplr(ul_laser)],'red','facealpha',0.3,'edgecolor','none')

title('Narrow-spiking single-unit activity around laser onset','fontweight','normal');
savefig(gcf,[folder filesep 'Mean Arch subject laser onset raster.fig']);

save([folder filesep 'Arch_subject_NS_onset.mat'],'psth_ctrl','psth_laser','T_vec');

%% Figures 3 and 5 and S7c (MU)

fields = {'perf','driven_FR','TS','RMS','tau'};
field_titles = {'Performance','Driven FR','Trial similarity','RMS difference','Optimal tau'};

SUorMU = {'single','mua'};

groups = {'single units','multi units'};

for g = 1:2

    tempstruct = mainset(strcmp({mainset.unit},SUorMU{g}));
    tempstruct = tempstruct([tempstruct.ctrl_perf] >= perf_thresh & [tempstruct.ctrl_d] >= d_thresh);

    for ff = 1:5 % measures

        perf_ctrl_clean = [tempstruct(strcmp({tempstruct.stim},'Clean')).(['ctrl_' fields{ff}])];
        perf_laser_clean = [tempstruct(strcmp({tempstruct.stim},'Clean')).(['laser_' fields{ff}])];

        perf_ctrl_masked = [tempstruct(strcmp({tempstruct.stim},'Masked')).(['ctrl_' fields{ff}])];
        perf_laser_masked = [tempstruct(strcmp({tempstruct.stim},'Masked')).(['laser_' fields{ff}])];

        [~,p_clean] = ttest(perf_ctrl_clean,perf_laser_clean)
        [~,p_masked] = ttest(perf_ctrl_masked,perf_laser_masked)
        [d_clean] = computeCohen_d(perf_ctrl_clean,perf_laser_clean,'paired')
        [d_masked] = computeCohen_d(perf_ctrl_masked,perf_laser_masked,'paired')

        figure('units','inches','position',[2 2 1.8 2.1]);
        subplot(1,2,1);
        plot((ones(numel(perf_ctrl_clean),2).*[1 2])',[perf_ctrl_clean;perf_laser_clean],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;
        set(gca,'xticklabels',{'Control','Laser'}); sigstar([1 2],p_clean); set(gca,'fontsize',8);
        set(gca,'xtick',[1 2]); title('Clean','Fontsize',10,'fontweight','normal'); % title(['Clean, p = ' num2str(p_clean)]);
        lim1 = get(gca,'ylim'); xlim([0.5 2.5]);

        if ff == 1, ytickformat('percentage'); end

        subplot(1,2,2);
        plot((ones(numel(perf_ctrl_masked),2).*[1 2])',[perf_ctrl_masked;perf_laser_masked],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;

        set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_masked);  set(gca,'fontsize',8);
        set(gca,'xtick',[1 2]); title('Masked','Fontsize',10,'fontweight','normal'); %title(['Masked, p = ' num2str(p_masked)]);
        lim2 = get(gca,'ylim'); xlim([0.5 2.5])

        if ff == 1, ytickformat('percentage'); end

        if ff == 1, lims = [50 104];
        else, lims = [lim1 lim2]; end

        subplot(1,2,1); ylim([min(lims) max(lims)]); plot([1 2],mean([perf_ctrl_clean;perf_laser_clean],2,'omitnan'),'-ok','linewidth',2,'markersize',4);
        subplot(1,2,2); ylim([min(lims) max(lims)]); plot([1 2],mean([perf_ctrl_masked;perf_laser_masked],2,'omitnan'),'-ok','linewidth',2,'markersize',4);

        sgtitle([field_titles{ff}],'fontsize',10);

        savefig(gcf,[folder filesep field_titles{ff} ', ' groups{g}]);
        [cc,p] = corrcoef([[tempstruct.ctrl_perf],[tempstruct.laser_perf]],[[tempstruct.(['ctrl_' fields{ff}])],[tempstruct.(['laser_' fields{ff}])]])

        pause;
        close;

    end
end

%% Plot optimal tau histogram - Figure 6a and 7b

taus = [[singleunits.ctrl_tau],[singleunits.laser_tau]]*1000;
bins = 0:10:250;
figure('unit','inches','position',[5 5 2.3 2]);
plot(5+bins(1:end-1),histcounts(taus,bins),'k','linewidth',1);
hold on;
plot(median(taus)*ones(1,2),[0;25],'--k');
Q = quantile(taus,[0.25 0.75]);
patch([Q(1) Q(2) Q(2) Q(1)],[0 0 25 25],[0.5 0.5 0.5],'edgecolor','none','facealpha',0.6); xlim([0 250]); box off;
xlabel('Optimal \tau (ms)'); ylabel('Count')
savefig(gcf,[folder filesep 'Optimal taus, single units.fig']);

%% VR-based performances for remaining single units - Figure 6b-d and Table 1
Tau = [0.0010    0.0020    0.0040    0.0080    0.0160    0.0320    0.0640    0.1280    0.2560]*1000;

p_clean = [];
p_masked = [];
d_clean = [];
d_masked = [];

figure('units','inches','position',[2 2 2.3 2.6]);

for ff = 1:length(Tau)
    perf_ctrl_clean = [tempstruct(strcmp({tempstruct.stim},'Clean')).(['ctrl_perf_' num2str(Tau(ff))])];
    perf_laser_clean = [tempstruct(strcmp({tempstruct.stim},'Clean')).(['laser_perf_' num2str(Tau(ff))])];

    perf_ctrl_masked = [tempstruct(strcmp({tempstruct.stim},'Masked')).(['ctrl_perf_' num2str(Tau(ff))])];
    perf_laser_masked = [tempstruct(strcmp({tempstruct.stim},'Masked')).(['laser_perf_' num2str(Tau(ff))])];

    [~,p_clean(ff)] = ttest(perf_ctrl_clean,perf_laser_clean)
    [~,p_masked(ff)] = ttest(perf_ctrl_masked,perf_laser_masked)
    d_clean(ff) = computeCohen_d(perf_ctrl_clean,perf_laser_clean, 'paired')
    d_masked(ff) = computeCohen_d(perf_ctrl_masked,perf_laser_masked, 'paired')

    subplot(1,2,1);
    plot((ones(numel(perf_ctrl_clean),2).*[1 2])',[perf_ctrl_clean;perf_laser_clean],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;
    set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_clean(ff)); set(gca,'fontsize',8);
    set(gca,'xtick',[1 2]); title('Clean','Fontsize',10); % title(['Clean, p = ' num2str(p_clean)]);
    xlim([0.5 2.5]);
    ytickformat('percentage');

    subplot(1,2,2);
    plot((ones(numel(perf_ctrl_masked),2).*[1 2])',[perf_ctrl_masked;perf_laser_masked],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;

    set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_masked(ff));  set(gca,'fontsize',8);
    set(gca,'xtick',[1 2]); title('Masked','Fontsize',10); %title(['Masked, p = ' num2str(p_masked)]);
    xlim([0.5 2.5]);
    ytickformat('percentage');

    lims = [50 104];

    subplot(1,2,1); ylim([min(lims) max(lims)]); plot([1 2],mean([perf_ctrl_clean;perf_laser_clean],2,'omitnan'),'-ok','linewidth',2,'markersize',4);
    subplot(1,2,2); ylim([min(lims) max(lims)]); plot([1 2],mean([perf_ctrl_masked;perf_laser_masked],2,'omitnan'),'-ok','linewidth',2,'markersize',4);

    bigtitle = ['VR-based performance @ ' num2str(Tau(ff)) 'ms'];

    sgtitle(bigtitle,'fontsize',10); disp(Tau(ff));
    %pause;

    if ismember(Tau(ff),[8 32 256])
        %savefig(gcf,[folder filesep bigtitle]);
    end
    clf;
end

VR_results = table(Tau',p_clean',p_masked',d_clean',d_masked','VariableNames',{'Tau','P_clean','P_masked','d_clean','d_masked'})
save([folder filesep 'VR_results.mat'],'VR_results');
%% Different performance metrics - Figure 4a-c

% implement logic for control and emergent hotspots
d1 = find( [perf_comps.ctrl_d] >= d_thresh & [perf_comps.ctrl_perf_SPIKE] >= perf_thresh);  % control

mainset_comp = perf_comps(d1);

singleunits_comp = mainset_comp(strcmp({mainset_comp.unit},'single'));

% different performances
real_comp = mainset_comp(structinds);

% load('mainset_diffperfs.mat');

perf_fields = {'ISI','RISPIKE','spkct'};
perf_titles = {'ISI-distance','RI-SPIKE-distance','Spike count distance'};
SUorMU = {'single','mua'};

for g = 1

    real_comp = mainset_comp(strcmp({mainset_comp.unit},SUorMU{g}));

    temp = struct2cell(real_comp);
    fields = fieldnames(real_comp);

    inds_clean = find(strcmp({real_comp.stim},'Clean'));
    inds_masked = find(strcmp({real_comp.stim},'Masked'));

    figure('units','inches','position',[2 2 2 2.3]);
    for ff = 1:3

        cf_ind = find(strcmp(fields,['ctrl_perf_' perf_fields{ff}]));
        lf_ind = find(strcmp(fields,['laser_perf_' perf_fields{ff}]));

        temp_ctrl_clean = vertcat(temp{cf_ind,:,inds_clean});
        temp_laser_clean = vertcat(temp{lf_ind,:,inds_clean});

        temp_ctrl_masked = vertcat(temp{cf_ind,:,inds_masked});
        temp_laser_masked = vertcat(temp{lf_ind,:,inds_masked});

        [~,p_clean(ff)] = ttest(temp_ctrl_clean,temp_laser_clean)
        [~,p_masked(ff)] = ttest(temp_ctrl_masked,temp_laser_masked)
        d_clean(ff) = computeCohen_d(temp_ctrl_clean,temp_laser_clean, 'paired')
        d_masked(ff) = computeCohen_d(temp_ctrl_masked,temp_laser_masked, 'paired')

        subplot(1,2,1);
        plot((ones(numel(temp_ctrl_clean),2).*[1 2])',[temp_ctrl_clean,temp_laser_clean]','-o','color',[0.86 0.86 0.86],'markersize',4);
        hold on;
        set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_clean(ff)); set(gca,'Fontsize',8);
        set(gca,'xtick',[1 2]); title('Clean','Fontsize',10,'Fontweight','normal'); % title(['Clean, p = ' num2str(p_clean)]);
        lim1 = get(gca,'ylim'); xlim([0.5 2.5]); ytickformat('percentage');

        subplot(1,2,2);
        plot((ones(numel(temp_ctrl_masked),2).*[1 2])',[temp_ctrl_masked,temp_laser_masked]','-o','color',[0.86 0.86 0.86],'markersize',4);
        hold on;

        set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_masked(ff));  set(gca,'Fontsize',8);
        set(gca,'xtick',[1 2]); title('Masked','Fontsize',10,'Fontweight','normal'); %title(['Masked, p = ' num2str(p_masked)]);
        lim2 = get(gca,'ylim'); xlim([0.5 2.5]); ytickformat('percentage');

        lims = [lim1 lim2];
        lims = [46 104];

        subplot(1,2,1); ylim([min(lims) max(lims)]); plot([1 2],mean([temp_ctrl_clean,temp_laser_clean]),'-ok','linewidth',2,'markersize',4);
        subplot(1,2,2); ylim([min(lims) max(lims)]); plot([1 2],mean([temp_ctrl_masked,temp_laser_masked]),'-ok','linewidth',2,'markersize',4);

        title1 = [perf_titles{ff} ' performance, ' groups{1+g}];
        sgtitle(title1);

        savefig(gcf,[folder filesep perf_titles{ff} ', ' groups{1+g}])

        pause;
        clf;
    end
    close;
end

%% Summary plot - Figure 4d

real_comp = mainset_comp(strcmp({mainset_comp.unit},'single')); % SU only

fields = fieldnames(real_comp);
temp = struct2cell(real_comp);

inds_clean = find(strcmp({real_comp.stim},'Clean'));
inds_masked = find(strcmp({real_comp.stim},'Masked'));

cf_SPK = find(strcmp(fields,'ctrl_perf_SPIKE'));
lf_SPK = find(strcmp(fields,'laser_perf_SPIKE'));

cf_RISPK = find(strcmp(fields,'ctrl_perf_RISPIKE'));
lf_RISPK = find(strcmp(fields,'laser_perf_RISPIKE'));

cf_ISI = find(strcmp(fields,'ctrl_perf_ISI'));
lf_ISI = find(strcmp(fields,'laser_perf_ISI'));

cf_spkct = find(strcmp(fields,'ctrl_perf_spkct'));
lf_spkct = find(strcmp(fields,'laser_perf_spkct'));

dists = {'SPK','ISI','RISPK','spkct'};

figure('units','inches','position',[3 3 4.8 2.4])
inputs_ctrl_clean = []; inputs_ctrl_masked = [];
inputs_laser_clean = []; inputs_laser_masked = [];
for d = 1:4
    x1 = [-0.2 0.2] + d;
    
    % clean
    subplot(1,2,1);
    temp_ctrl_clean = vertcat(temp{eval(['cf_' dists{d}]),:,inds_clean});
    temp_laser_clean = vertcat(temp{eval(['lf_' dists{d}]),:,inds_clean});

    hold on;
    scatter(ones(size(temp_ctrl_clean)) * x1(1),temp_ctrl_clean,36,'Marker','o','markeredgecolor',[0.86 0.86 0.86]);
    scatter(ones(size(temp_laser_clean)) * x1(2),temp_laser_clean,49,'Marker','s','markeredgecolor',[0.86 0.86 0.86]);
    plot(x1,[temp_ctrl_clean';temp_laser_clean'],'color',[0.86 0.86 0.86]);
    
    % mean
    scatter(x1(1),mean(temp_ctrl_clean),36,'k','Marker','o','linewidth',2);
    scatter(x1(2),mean(temp_laser_clean),49,'k','Marker','s','linewidth',2);
    plot(x1,mean([temp_ctrl_clean,temp_laser_clean]),'-k','linewidth',2);

    xlim([0.5 4.5]);
    ylim([46 104]); ytickformat('percentage'); set(gca,'xticklabel',{'SPIKE','ISI','RI-SPIKE','Spk ct'},'xtick',1:4);
    title('Clean','Fontweight','normal','Fontsize',10); ylabel('Performance');
    set(gca,'Fontsize',8);

    subplot(1,2,2);
    temp_ctrl_masked = vertcat(temp{eval(['cf_' dists{d}]),:,inds_masked});
    temp_laser_masked = vertcat(temp{eval(['lf_' dists{d}]),:,inds_masked});

    hold on;
    scatter(ones(size(temp_ctrl_masked)) * x1(1),temp_ctrl_masked,36,'Marker','o','markeredgecolor',[0.86 0.86 0.86]);
    scatter(ones(size(temp_laser_masked)) * x1(2),temp_laser_masked,49,'Marker','s','markeredgecolor',[0.86 0.86 0.86]);
    plot(x1,[temp_ctrl_masked';temp_laser_masked'],'color',[0.86 0.86 0.86]);
    
    % mean
    scatter(x1(1),mean(temp_ctrl_masked),36,'k','Marker','o','linewidth',2);
    scatter(x1(2),mean(temp_laser_masked),49,'k','Marker','s','linewidth',2);
    plot(x1,mean([temp_ctrl_masked,temp_laser_masked]),'-k','linewidth',2);

    xlim([0.5 4.5]);
    ylim([46 104]);  ytickformat('percentage'); set(gca,'xticklabel',{'SPIKE','ISI','RI-SPIKE','Spk ct'},'xtick',1:4);
    title('Masked','Fontweight','normal','Fontsize',10);
    set(gca,'Fontsize',8);

    inputs_ctrl_clean = cat(2,inputs_ctrl_clean,temp_ctrl_clean);
    inputs_ctrl_masked = cat(2,inputs_ctrl_masked,temp_ctrl_masked);

    inputs_laser_clean = cat(2,inputs_laser_clean,temp_laser_clean);
    inputs_laser_masked = cat(2,inputs_laser_masked,temp_laser_masked);

end

% copy sigstars
for d = 1:4
    x1 = [-0.2 0.2] + d;

    subplot(1,2,1);
    [~,p] = ttest(inputs_ctrl_clean(:,d),inputs_laser_clean(:,d));
    sigstar(x1,p)

    subplot(1,2,2);
    [~,p] = ttest(inputs_ctrl_masked(:,d),inputs_laser_masked(:,d));
    sigstar(x1,p)
end

savefig(gcf,[folder filesep 'Performance bar plots, single units.fig']);

%% Upper envelope - Figure 2f

tempstruct = mainset(strcmp({mainset.unit},'single')); % SU only

for ns = 1:length(exps)
    good_chans{ns} = unique([tempstruct(strcmp({tempstruct.subject},subjects{ns})).ch]);
end

plotPopulationEnv_V2(exps,good_chans,whichTMR,whichpower,folder)

%% Figure S2c - CSD sink onset

[mean_onset] = plotOnsetDepth(exps,whichpower,whichTMR);

%% optimal tau histogram, control AND laser - Figure 6

tempstruct = mainset(strcmp({mainset.unit},'single')); % SU only

taus = [[tempstruct.ctrl_tau],[tempstruct.laser_tau]]*1000;
bins = 0:10:250;
figure('unit','inches','position',[5 5 2.3 2]);
plot(5+bins(1:end-1),histcounts(taus,bins),'k','linewidth',1);
hold on;
plot(median(taus)*ones(1,2),[0;25],'--k');
Q = quantile(taus,[0.25 0.75]);
patch([Q(1) Q(2) Q(2) Q(1)],[0 0 25 25],[0.5 0.5 0.5],'edgecolor','none','facealpha',0.6); xlim([0 250]); box off;
xlabel('Optimal \tau (ms)'); ylabel('Count')
savefig(gcf,[folder filesep 'Optimal taus, single units.fig']);


%% Fig S7 for layers and single unit NS vs RS

% run anova comparisons for SU/MU and RS/NS

meas = table({'Ctrl','Laser'}','VariableNames',{'Cond'});
meas.Cond = categorical(meas.Cond);

tempstruct = mainset(strcmp({mainset.unit},'single')); % SU only

clean_results = tempstruct(strcmp({tempstruct.stim},'Clean'));
masked_results = tempstruct(strcmp({tempstruct.stim},'Masked'));

perf_ctrl_clean = [clean_results.ctrl_perf]';
perf_laser_clean = [clean_results.laser_perf]';

perf_ctrl_masked = [masked_results.ctrl_perf]';
perf_laser_masked = [masked_results.laser_perf]';

%%%%%%%%%% RS vs NS

clean_units = strcmp({clean_results.type},'narrow_spiking')';
masked_units = strcmp({masked_results.type},'narrow_spiking')';

laser_inds = [zeros(numel(perf_laser_clean),1);ones(numel(perf_laser_clean),1)];

[~,inds] = sort(clean_units);
X = [perf_ctrl_clean(inds);perf_laser_clean(inds)];
group = [laser_inds,repmat(clean_units(inds),2,1)];
stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','unit'});

laser_inds = [zeros(numel(perf_ctrl_masked),1);ones(numel(perf_ctrl_masked),1)];

[~,inds] = sort(masked_units);
X = [perf_ctrl_masked(inds);perf_laser_masked(inds)];
group = [laser_inds,repmat(masked_units(inds),2,1)];
stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','unit'});

% clean
tab = table(perf_ctrl_clean,perf_laser_clean,clean_units,'VariableNames',{'Ctrl','Laser','NarrowSpk'});
tab = convertvars(tab,'NarrowSpk','categorical');
rm = fitrm(tab,'Ctrl-Laser~NarrowSpk','WithinDesign',meas);
tbl = ranova(rm,'WithinModel','Cond')
within_layers = multcompare(rm,'NarrowSpk','by','Cond'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
within_layers = multcompare(rm,'Cond','by','NarrowSpk'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
etaSq_layer = tbl{2,4}*tbl{2,2} / (tbl{2,4}*tbl{2,2} + tbl{3,2})
etaSq_int = tbl{5,4}*tbl{5,2} / (tbl{5,4}*tbl{5,2} + tbl{6,2})


% masked
tab = table(perf_ctrl_masked,perf_laser_masked,masked_units,'VariableNames',{'Ctrl','Laser','NarrowSpk'});
tab = convertvars(tab,'NarrowSpk','categorical');
rm = fitrm(tab,'Ctrl-Laser~NarrowSpk','WithinDesign',meas);
tbl = ranova(rm,'WithinModel','Cond')
within_layers = multcompare(rm,'NarrowSpk','by','Cond'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
within_layers = multcompare(rm,'Cond','by','NarrowSpk'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
etaSq_layer = tbl{2,4}*tbl{2,2} / (tbl{2,4}*tbl{2,2} + tbl{3,2})
etaSq_int = tbl{5,4}*tbl{5,2} / (tbl{5,4}*tbl{5,2} + tbl{6,2})

%%%%%%%%%%%%%%%%%%%%%%%% Layers

layers_clean = zeros(size(perf_ctrl_clean));
layers_masked = zeros(size(perf_ctrl_masked));

for n = 2:4
    layers_clean(strcmp({clean_results.layer},layer_names{n})) = n-1;%layer_names{n};
    layers_masked(strcmp({masked_results.layer},layer_names{n})) = n-1;%layer_names{n};
end

laser_inds = [zeros(numel(perf_ctrl_clean),1);ones(numel(perf_ctrl_clean),1)];

[~,inds] = sort(layers_clean);
X = [perf_ctrl_clean(inds);perf_laser_clean(inds)];
group = [laser_inds,repmat(layers_clean(inds),2,1)];
stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','layer'});

laser_inds = [zeros(numel(perf_ctrl_masked),1);ones(numel(perf_ctrl_masked),1)];

[~,inds] = sort(layers_masked);
X = [perf_ctrl_masked(inds);perf_laser_masked(inds)];
group = [laser_inds,repmat(layers_masked(inds),2,1)];
stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','layer'});

% clean
tab = table(perf_ctrl_clean,perf_laser_clean,layers_clean,'VariableNames',{'Ctrl','Laser','Layer'});
tab = convertvars(tab,'Layer','categorical');
rm = fitrm(tab,'Ctrl-Laser~Layer','WithinDesign',meas);
tbl = ranova(rm,'WithinModel','Cond')
within_conds = multcompare(rm,'Layer','by','Cond'); [~,ind] = sort(within_conds{:,4}); within_conds = within_conds(ind,:); within_conds(1:size(within_conds,1)/2,:) = []
within_layers = multcompare(rm,'Cond','by','Layer'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
etaSq_layer = tbl{2,4}*tbl{2,2} / (tbl{2,4}*tbl{2,2} + tbl{3,2})
etaSq_int = tbl{5,4}*tbl{5,2} / (tbl{5,4}*tbl{5,2} + tbl{6,2})

% masked
tab = table(perf_ctrl_masked,perf_laser_masked,layers_masked,'VariableNames',{'Ctrl','Laser','Layer'});
tab = convertvars(tab,'Layer','categorical');
rm = fitrm(tab,'Ctrl-Laser~Layer','WithinDesign',meas);
tbl = ranova(rm,'WithinModel','Cond')
within_cond = multcompare(rm,'Layer','by','Cond'); [~,ind] = sort(within_cond{:,4}); within_cond = within_cond(ind,:); within_cond(1:size(within_cond,1)/2,:) = []
within_layers = multcompare(rm,'Cond','by','Layer'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
etaSq_layer = tbl{2,4}*tbl{2,2} / (tbl{2,4}*tbl{2,2} + tbl{3,2})
etaSq_int = tbl{5,4}*tbl{5,2} / (tbl{5,4}*tbl{5,2} + tbl{6,2})

%% S7 plots - Show performances for different group comparisons (SU vs MU, NS vs RS, Layers)

% for NS vs RS and layers, SINGLE UNITS ONLY

comps = {'unit','type','layer'};
booleans = {'single','narrow_spiking'};
type = {'SingleUnit','Spiking','Layer'};

% clean = mainset(strcmp({mainset.stim},'Clean'));
% masked = mainset(strcmp({mainset.stim},'Masked'));

% real = mainset(structinds == 1);

perf_fields = {'SPIKE','ISI','RISPIKE'};
groups = {{'SU','MU'},{'RS','NS'},layer_names};

for si = 1:length(structinds)
    if structinds(si) == 1
        mainset_comp(si).unit = 'single';
    else
        mainset_comp(si).unit = 'mua';
    end
end

% as barplots
figure('unit','inches','position',[4 4 3 2.5]);

% g: different types of groups
% 1: SU vs MU, 2: RS vs NS, 3: vs layers

for g = 2:3
    % focus on SU data only for the other two groupings
    if g > 1
        tempstruct = mainset_comp(structinds == 1); % SUs only
    else
        tempstruct = mainset_comp;
    end
    clean = tempstruct(strcmp({tempstruct.stim},'Clean'));
    masked = tempstruct(strcmp({tempstruct.stim},'Masked'));
    groupIDs = groups{g};

    % p: different spike distance metrics (see perf_fields)
    for p = 1

        input_clean_ctrl = cell(1,1);
        input_clean_laser = cell(1,1);
        input_masked_ctrl = cell(1,1);
        input_masked_laser = cell(1,1);

        % group mainset results by type
        if g < 3
            % clean
            inds_clean = strcmp({clean.(comps{g})},booleans{g});

            % flip logic for NS/RS comparisons so that NS data is x =2 in
            % barplots
            if g == 2
                inds_clean = ~inds_clean;
            end

            input_clean_ctrl{1} = [clean(inds_clean).(['ctrl_perf_' perf_fields{p}])];
            input_clean_laser{1} = [clean(inds_clean).(['laser_perf_' perf_fields{p}])];
            input_clean_ctrl{2} = [clean(~inds_clean).(['ctrl_perf_' perf_fields{p}])];
            input_clean_laser{2} = [clean(~inds_clean).(['laser_perf_' perf_fields{p}])];

            % masked
            inds_masked = strcmp({masked.(comps{g})},booleans{g});
            if g == 2
                inds_masked = ~inds_masked;
            end

            input_masked_ctrl{1} = [masked(inds_masked).(['ctrl_perf_' perf_fields{p}])];
            input_masked_laser{1} = [masked(inds_masked).(['laser_perf_' perf_fields{p}])];
            input_masked_ctrl{2} = [masked(~inds_masked).(['ctrl_perf_' perf_fields{p}])];
            input_masked_laser{2} = [masked(~inds_masked).(['laser_perf_' perf_fields{p}])];

            % RS,SU == 1 and NS,MU == 2 in groups
            inds_clean = mod(inds_clean+1,2) + 1;
            inds_masked = mod(inds_masked+1,2) + 1;
        else
            % compile clean and masked results per layer
            for l = 1:5
                input_clean_ctrl{l} = [clean([clean.layer] == l).(['ctrl_perf_' perf_fields{p}])];
                input_clean_laser{l} = [clean([clean.layer] == l).(['laser_perf_' perf_fields{p}])];

                input_masked_ctrl{l} = [masked([masked.layer] == l).(['ctrl_perf_' perf_fields{p}])];
                input_masked_laser{l} = [masked([masked.layer] == l).(['laser_perf_' perf_fields{p}])];
            end
            inds_clean = [clean.layer];
            inds_masked = [masked.layer];
            set(gcf,'position',[4 4 4.2 3]);
        end
        x1 = 1:numel(input_clean_ctrl);

        % clean
        subplot(1,2,1); hold on;
        bar(x1 - 0.2,cellfun(@mean,input_clean_ctrl),0.4,'facecolor','none','linewidth',2); % control
        errorbar(x1 - 0.2,cellfun(@mean,input_clean_ctrl),cellfun(@(x) std(x)/sqrt(numel(x)),input_clean_ctrl),'k','linestyle','none')
        bar(x1 + 0.2,cellfun(@mean,input_clean_laser),0.4,'facecolor','k','linewidth',2); % laser
        errorbar(x1 + 0.2,cellfun(@mean,input_clean_laser),cellfun(@(x) std(x)/sqrt(numel(x)),input_clean_laser),'k','linestyle','none')
        if g < 3, xlim([0.5 x1(end)+0.5]); else, xlim([1.5 x1(end)-0.5]); end
        ylim([50 90]); set(gca,'fontsize',8,'xtick',x1,'xticklabel',groups{g}); ytickformat('percentage');
        title('Clean','fontweight','normal','fontsize',10);

        % masked
        subplot(1,2,2); hold on;
        h1 = bar(x1 - 0.2,cellfun(@mean,input_masked_ctrl),0.4,'facecolor','none','linewidth',2); % control
        errorbar(x1 - 0.2,cellfun(@mean,input_masked_ctrl),cellfun(@(x) std(x)/sqrt(numel(x)),input_masked_ctrl),'k','linestyle','none')
        h2 = bar(x1 + 0.2,cellfun(@mean,input_masked_laser),0.4,'facecolor','k','linewidth',2); % laser
        errorbar(x1 + 0.2,cellfun(@mean,input_masked_laser),cellfun(@(x) std(x)/sqrt(numel(x)),input_masked_laser),'k','linestyle','none')
        if g < 3, xlim([0.5 x1(end)+0.5]); else, xlim([1.5 x1(end)-0.5]); end
        ylim([50 90]); set(gca,'fontsize',8,'xtick',x1,'xticklabel',groups{g}); %,'ytick',[]);
        title('Masked','fontweight','normal','fontsize',10);

        sgtitle([perf_fields{p} ' performance']);

        % clean ANOVA
        meas = table({'Ctrl','Laser'}','VariableNames',{'Cond'});
        meas.Cond = categorical(meas.Cond);
        tab = table([clean.(['ctrl_perf_' perf_fields{p}])]',[clean.(['laser_perf_' perf_fields{p}])]',groupIDs(inds_clean)','VariableNames',{'Ctrl','Laser',type{g}});
        tab = convertvars(tab,type{g},'categorical');
        rm = fitrm(tab,['Ctrl-Laser~' type{g}],'WithinDesign',meas);
        tbl = ranova(rm,'WithinModel','Cond')

        %%%%% For eta-squared calculations
        temp_input_c = [clean.(['ctrl_perf_' perf_fields{p}])]';
        temp_input_l = [clean.(['laser_perf_' perf_fields{p}])]';

        laser_inds = [zeros(numel(temp_input_c),1);ones(numel(temp_input_c),1)];
        [~,inds] = sort(inds_clean); % sort by increasing layer
        X = [temp_input_c(inds);temp_input_l(inds)];
        group = [laser_inds,repmat(inds_clean(inds)',2,1)];
        stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','layer'})


        anova_tbl = formattedDisplayText(tbl);

        if g == 3, x1([1 end]) = []; end

        if any(tbl{2:3:end,5} < 0.05)
            within_conds = multcompare(rm,type{g},'by','Cond'); [~,ind] = sort(within_conds{:,4}); within_conds = within_conds(ind,:); within_conds(1:size(within_conds,1)/2,:) = [];
            within_groups = multcompare(rm,'Cond','by',type{g}); [~,ind] = sort(within_groups{:,4}); within_groups = within_groups(ind,:); within_groups(1:size(within_groups,1)/2,:) = [];

            subplot(1,2,1);
            % SU / RS should be first group
            if ~strcmp(groups{g}{1},cellstr(within_groups{1,1}))
                within_groups = flipud(within_groups)
            else
                within_groups
            end
            if g == 3, within_groups = sortrows(within_groups); end
            pvals = within_groups{:,6}; grps = num2cell([x1-0.2;x1+0.2],1);
            sigstar(grps,pvals);

            % control should be first group
            within_conds = sortrows(within_conds,[1 3])
            pvals = within_conds{:,6}; grps = num2cell([x1-0.2;x1+0.2]',1);
            if g == 3, grps = num2cell([nchoosek(x1-0.2,2);nchoosek(x1+0.2,2)],2); end
            sigstar(grps,pvals);

        else
            within_conds = 0; within_groups = 0;
        end
        save([folder filesep perf_fields{p} '_perf_' type{g} '_clean_stats.mat'],'anova_tbl','within_conds','within_groups');

        % masked ANOVA
        meas = table({'Ctrl','Laser'}','VariableNames',{'Cond'});
        meas.Cond = categorical(meas.Cond);
        tab = table([masked.(['ctrl_perf_' perf_fields{p}])]',[masked.(['laser_perf_' perf_fields{p}])]',groupIDs(inds_masked)','VariableNames',{'Ctrl','Laser',type{g}});
        tab = convertvars(tab,type{g},'categorical');
        rm = fitrm(tab,['Ctrl-Laser~' type{g}],'WithinDesign',meas);
        tbl = ranova(rm,'WithinModel','Cond')

        %%%%% For eta-squared calculations
        temp_input_c = [masked.(['ctrl_perf_' perf_fields{p}])]';
        temp_input_l = [masked.(['laser_perf_' perf_fields{p}])]';

        laser_inds = [zeros(numel(temp_input_c),1);ones(numel(temp_input_c),1)];
        [~,inds] = sort(inds_masked); % sort by increasing layer
        X = [temp_input_c(inds);temp_input_l(inds)];
        group = [laser_inds,repmat(inds_masked(inds)',2,1)];
        stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','layer'})


        anova_tbl = formattedDisplayText(tbl);

        if any(tbl{2:3:end,5} < 0.05)
            within_conds = multcompare(rm,type{g},'by','Cond'); [~,ind] = sort(within_conds{:,4}); within_conds = within_conds(ind,:); within_conds(1:size(within_conds,1)/2,:) = []
            within_groups = multcompare(rm,'Cond','by',type{g}); [~,ind] = sort(within_groups{:,4}); within_groups = within_groups(ind,:); within_groups(1:size(within_groups,1)/2,:) = []

            subplot(1,2,2);
            % SU / RS should be first group
            if ~strcmp(groups{g}{1},cellstr(within_groups{1,1}))
                within_groups = flipud(within_groups)
            end
            if g == 3, within_groups = sortrows(within_groups); end
            pvals = within_groups{:,6}; grps = num2cell([x1-0.2;x1+0.2],1);
            sigstar(grps,pvals);

            within_conds = sortrows(within_conds,[1 3]);
            % control should be first group
            pvals = within_conds{:,6}; grps = num2cell([x1-0.2;x1+0.2]',1);
            if g == 3, grps = num2cell([nchoosek(x1-0.2,2);nchoosek(x1+0.2,2)],2); end
            sigstar(grps,pvals);
        else
            within_conds = 0; within_groups = 0;
        end
        save([folder filesep  perf_fields{p} '_perf_' type{g} '_masked_stats.mat'],'anova_tbl','within_conds','within_groups');

        legend([h1 h2],'Control','Laser','location','best');
        saveas(gcf,[folder filesep  perf_fields{p} ' performance, ' type{g} ' comparisons.png']);
        savefig(gcf,[folder filesep  perf_fields{p} ' performance, ' type{g} ' comparisons']);
        pause;
        clf;

    end
end

%% S6 - Laser effect vs layer for RS and NS in single units only

layer_names = {'L1','L2/3','L4','L5','L6'};
layer_edges = [0.5 3.5 6.5 8.5 11.5 15.5];

chanMap = [9 8 6 11 7 10 4 12;
    13 1 15 3 14 2 16 5;
    31 19 29 17 32 20 27 18;
    26 23 21 28 24 25 22 30];

n = 1; layer_results = struct;
for ns = 1:length(exps)

    Chans = SU{ns};

    for ch = Chans

        [shank,depth] = find(chanMap == ch);
        temp_dpth = 8 + depth - exps(ns).gran_layers(shank);
        layer_ind = discretize(temp_dpth,layer_edges);

        spk_width = calcSpkWidth(exps(ns).ctrl_spks.snips_clean{ch},exps(ns).ctrl_spks.snips_masked{ch});
        if spk_width < 0.5
            celltype = 'NS';
            % continue
        else
            celltype = 'RS';
        end

        % clean
        ctrl_driven_FR = []; laser_driven_FR = []; ctrl_spon_FR = []; laser_spon_FR = []; ctrl_onset_FR = []; laser_onset_FR = [];
        for i = 1:4
            spks_ctrl = squeeze(exps(ns).ctrl_spks.Spks_clean{ch}(:,i,:));
            spks_laser = squeeze(exps(ns).laser_spks.Spks_clean{ch}(:,i,:));

            ctrl_driven_FR = cat(1,ctrl_driven_FR,calcFRPeriod(spks_ctrl,'driven'));
            laser_driven_FR = cat(1,laser_driven_FR,calcFRPeriod(spks_laser,'driven'));
            ctrl_onset_FR = cat(1,ctrl_onset_FR,calcFRPeriod(spks_ctrl,'onset'));
            laser_onset_FR = cat(1,laser_onset_FR,calcFRPeriod(spks_laser,'onset'));
            ctrl_spon_FR = cat(1,ctrl_spon_FR,calcFRPeriod(spks_ctrl,'laser'));
            laser_spon_FR = cat(1,laser_spon_FR,calcFRPeriod(spks_laser,'laser'));
        end

        % masked
        for i = 1:16
            [tloc,mloc] = ind2sub([4 4],i);

            spks_ctrl = squeeze(exps(ns).ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:));
            spks_laser = squeeze(exps(ns).laser_spks.Spks_masked{ch}(:,tloc,mloc,:));

            ctrl_driven_FR = cat(1,ctrl_driven_FR,calcFRPeriod(spks_ctrl,'driven'));
            laser_driven_FR = cat(1,laser_driven_FR,calcFRPeriod(spks_laser,'driven'));
            ctrl_onset_FR = cat(1,ctrl_onset_FR,calcFRPeriod(spks_ctrl,'onset'));
            laser_onset_FR = cat(1,laser_onset_FR,calcFRPeriod(spks_laser,'onset'));
            ctrl_spon_FR = cat(1,ctrl_spon_FR,calcFRPeriod(spks_ctrl,'laser'));
            laser_spon_FR = cat(1,laser_spon_FR,calcFRPeriod(spks_laser,'laser'));
        end
        layer_results(n).ctrl_driven_FR = mean(ctrl_driven_FR);
        layer_results(n).laser_driven_FR = mean(laser_driven_FR);
        layer_results(n).ctrl_onset_FR = mean(ctrl_onset_FR);
        layer_results(n).laser_onset_FR = mean(laser_onset_FR);
        layer_results(n).ctrl_spon_FR = mean(ctrl_spon_FR);
        layer_results(n).laser_spon_FR = mean(laser_spon_FR);

        layer_results(n).layer = layer_ind;
        layer_results(n).type = celltype;
        layer_results(n).subject = subjects{ns};
        layer_results(n).ch = ch;
        n = n + 1;
    end
end


groups = {'RS','NS'};
measures = {'spon_FR','onset_FR','driven_FR'};
titles = {'Spontaneous FR','Onset FR','Evoked FR'};

meas = table({'Ctrl','Laser'}','VariableNames',{'Cond'});
meas.Cond = categorical(meas.Cond);

for g = 1:2 % 1:RS ,2:NS

    tempstruct = layer_results(strcmp({layer_results.type},groups{g}));

    for m = 1:2 % 1:spontaneous, 2:onset, 3:driven
        figure('unit','inches','position',[4 4 2.25 2]);

        ctrl_res = [tempstruct.(['ctrl_' measures{m}])];
        laser_res = [tempstruct.(['laser_' measures{m}])];
        inds = [tempstruct.layer];
        scatter(inds,ctrl_res,25,'+k');
        hold on;
        scatter(inds,laser_res,16,'or');
        mean_ctrl = []; mean_laser = [];

        % mean per layer
        for la = 1:5
            mean_ctrl(la) = mean(ctrl_res(inds == la));
            mean_laser(la) = mean(laser_res(inds == la));
        end

        line(1:5,mean_ctrl,'color','k','linewidth',2);
        line(1:5,mean_laser,'color','r','linewidth',2);

        set(gca,'fontsize',8,'xticklabels',layer_names,'xtick',1:5,'fontname','arial');
        xlim([0.5 5.5]); ylim([0 60]);
        legend('Control','Laser'); ylabel('FR (Hz)'); xlabel('Layer');
        title([titles{m} ', ' groups{g} ' units'],'fontweight','normal');

        saveas(gcf,[folder filesep titles{m} ' vs layer, ' groups{g} ' units.png']);
        savefig(gcf,[folder filesep titles{m} ' vs layer, ' groups{g} ' units']);

        % remove L1 data
        ctrl_res(inds == 1) = [];
        laser_res(inds == 1) = [];
        inds(inds == 1) = [];

        layers = layer_names(inds);

        tab = table(ctrl_res',laser_res',layers','VariableNames',{'Ctrl','Laser','Layer'});
        tab = convertvars(tab,'Layer','categorical');
        rm = fitrm(tab,'Ctrl-Laser~Layer','WithinDesign',meas);
        tbl = ranova(rm,'WithinModel','Cond')

        laser_inds = [zeros(numel(ctrl_res),1);ones(numel(laser_res),1)];

        [~,is] = sort(inds);
        X = [ctrl_res(is)';laser_res(is)'];
        group = [laser_inds,repmat(inds(is)',2,1)];
        stats = mes2way(X,group,'partialeta2','isDep',[1 0],'fName',{'laser','unit'});
        %stats.partialeta2

        %         if any(tbl{2:3:end,5} < 0.05)
        within_conds = multcompare(rm,'Layer','by','Cond'); [~,ind] = sort(within_conds{:,4}); within_conds = within_conds(ind,:); within_conds(1:size(within_conds,1)/2,:) = []
        within_layers = multcompare(rm,'Cond','by','Layer'); [~,ind] = sort(within_layers{:,4}); within_layers = within_layers(ind,:); within_layers(1:size(within_layers,1)/2,:) = []
        %         end
        % save([folder filesep titles{m} '_layer_' groups{g} '_units'],'tbl','within_layers','within_conds');
        pause;
    end
end

