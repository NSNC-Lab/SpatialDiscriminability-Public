
clear; close all; clc

addpath('stimuli');
addpath(genpath('subfunctions'));
addpath('sigstar-master');
addpath(genpath('sortingQuality-master'));

analysis_storage = 'Control-data/analysis'; %uigetdir(cd,'Select where to save processed spikes and results.');   % 
data_storage = 'Control-data/processed-Kilosort'; %uigetdir(cd,'Select folder with stored data tanks.');     

addpath(analysis_storage);

perf_thresh = 70; d_thresh = 1;
whichpower = 10; whichTMR = 0;

subjects = {'150R-230201','155L-230126','151L-230206','152R-230206','153L-230206'}; %10mW, 0dB TMR

exps = struct;

perf_methods = {'_SPIKE'};

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

    load(fullfile(subfolder(expind).folder,subfolder(expind).name,'unitInfo.mat'));
    exps(ns).unitInfo = unitInfo;
end
fs = 24414.0625; clc

folder = sprintf('Control group %imW %idB TMR',whichpower,whichTMR);
if ~isfolder(folder); mkdir(folder); end

for ns = 1:length(exps)
    temp = exps(ns).unitInfo{:,1}; chs{ns} = temp';
    isSU = logical(exps(ns).unitInfo{:,5}); SU{ns} = temp(isSU);
    isNS = logical(exps(ns).unitInfo{:,6}); NS{ns} = temp(isNS);
end

%% Find hotspots

weak_chs = {[13 17],[29],[],[27],[]}; % 150R date 1/25, 155L date 1/26

locs = [90 45 0 -90];
results = struct; n = 1; unitnum = 1;

for ns = 1:length(exps)

    Chans = findResponsiveChans(exps(ns).ctrl_spks);
    Chans = setdiff(Chans,weak_chs{ns});

    for ch = Chans

        if ismember(ch,SU{ns}), mua = 'single'; else, mua = 'mua'; end
        if ismember(ch,NS{ns}), celltype = 'narrow_spiking'; else, celltype = 'regular'; end

        % check if performance is from low # trials
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end

        inds = [];

        % clean
        inds_ctrl = findHotspots(exps(ns).ctrl_perf_SPIKE.Max{ch},exps(ns).ctrl_perf_SPIKE.d{ch},...
            exps(ns).ctrl_spks.Spks_clean{ch},perf_thresh,d_thresh);
        inds_laser = findHotspots(exps(ns).laser_perf_SPIKE.Max{ch},exps(ns).laser_perf_SPIKE.d{ch},...
            exps(ns).laser_spks.Spks_clean{ch},perf_thresh,d_thresh);
        inds = union(inds_ctrl,inds_laser);
        if isrow(inds), inds = inds'; end

        if ~isempty(inds)
            for ii = 1:length(inds)

                spks_ctrl = squeeze(exps(ns).ctrl_spks.Spks_clean{ch}(:,inds(ii),:));
                spks_laser = squeeze(exps(ns).laser_spks.Spks_clean{ch}(:,inds(ii),:));

                temp_ctrl = []; temp_laser = [];

                % store spk trains in 2 cells, 1 per target
                for t = 1:2
                    trials = find(~cellfun(@isempty,spks_ctrl(:,t)));
                    nonemp_ctrl{t} = spks_ctrl(trials,t);

                    trials = find(~cellfun(@isempty,spks_laser(:,t)));
                    nonemp_laser{t} = spks_laser(trials,t);
                end
                minlen_ctrl = min(cellfun(@length,nonemp_ctrl));
                minlen_laser = min(cellfun(@length,nonemp_laser));

                for t = 1:2
                    nonemp_ctrl{t} = nonemp_ctrl{t}(1:minlen_ctrl)';
                    nonemp_laser{t} = nonemp_laser{t}(1:minlen_laser)';
                end

                % concatenate spk_trains and convert into a psth for
                % rcorr and sparseness functions
                for t = 1:minlen_ctrl
                    train1 = nonemp_ctrl{1}{t};
                    train2 = nonemp_ctrl{2}{t};
                    temp_ctrl(t,:) = histcounts([train1(train1 >=0 & train1 <= 3);...
                        train2(train2 >=0 & train2 <= 3)+3]',0:1/fs:6);
                end

                for t = 1:minlen_laser
                    train1 = nonemp_laser{1}{t};
                    train2 = nonemp_laser{2}{t};
                    temp_laser(t,:) = histcounts([train1(train1 >=0 & train1 <= 3);...
                        train2(train2 >=0 & train2 <= 3)+3]',0:1/fs:6);
                end

                results(n).ctrl_perf = exps(ns).ctrl_perf_SPIKE.Max{ch}(inds(ii));
                results(n).laser_perf = exps(ns).laser_perf_SPIKE.Max{ch}(inds(ii));
                results(n).ctrl_d = exps(ns).ctrl_perf_SPIKE.d{ch}(inds(ii));
                results(n).laser_d = exps(ns).laser_perf_SPIKE.d{ch}(inds(ii));
                results(n).ctrl_driven_FR = calcFRPeriod(spks_ctrl,'driven');
                results(n).laser_driven_FR = calcFRPeriod(spks_laser,'driven');
                results(n).ctrl_spon_FR = calcFRPeriod(spks_ctrl,'laser');
                results(n).laser_spon_FR = calcFRPeriod(spks_laser,'laser');
                results(n).unit = mua;
                results(n).stim = 'Clean';
                results(n).ch = ch;
                results(n).tloc = locs(inds(ii));
                results(n).mloc = nan;
                results(n).subject = subjects{ns};
                results(n).num = unitnum;
                results(n).celltype = celltype;
                n = n + 1;
            end
        end

        inds = [];

        % masked
        inds_ctrl = findHotspots(exps(ns).ctrl_perf_SPIKE.max_masked{ch},exps(ns).ctrl_perf_SPIKE.d_masked{ch},...
            exps(ns).ctrl_spks.Spks_masked{ch},perf_thresh,d_thresh);
        inds_laser = findHotspots(exps(ns).laser_perf_SPIKE.max_masked{ch},exps(ns).laser_perf_SPIKE.d_masked{ch},...
            exps(ns).laser_spks.Spks_masked{ch},perf_thresh,d_thresh);
        inds = union(inds_ctrl,inds_laser);
        if isrow(inds), inds = inds'; end

        if ~isempty(inds)
            for ii = 1:length(inds)
                [mloc,tloc] = ind2sub([4 4],inds(ii));
                mloc = 5 - mloc;

                spks_ctrl = squeeze(exps(ns).ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:));
                spks_laser = squeeze(exps(ns).laser_spks.Spks_masked{ch}(:,tloc,mloc,:));

                temp_ctrl = []; temp_laser = [];

                % store spk trains in 2 cells, 1 per target
                for t = 1:2
                    trials = find(~cellfun(@isempty,spks_ctrl(:,t)));
                    nonemp_ctrl{t} = spks_ctrl(trials,t);

                    trials = find(~cellfun(@isempty,spks_laser(:,t)));
                    nonemp_laser{t} = spks_laser(trials,t);
                end
                minlen_ctrl = min(cellfun(@length,nonemp_ctrl));
                minlen_laser = min(cellfun(@length,nonemp_laser));

                for t = 1:2
                    nonemp_ctrl{t} = nonemp_ctrl{t}(1:minlen_ctrl)';
                    nonemp_laser{t} = nonemp_laser{t}(1:minlen_laser)';
                end

                % concatenate spk_trains and convert into a psth for
                % rcorr and sparseness functions
                for t = 1:minlen_ctrl
                    train1 = nonemp_ctrl{1}{t};
                    train2 = nonemp_ctrl{2}{t};
                    temp_ctrl(t,:) = histcounts([train1(train1 >=0 & train1 <= 3);...
                        train2(train2 >=0 & train2 <= 3)+3]',0:1/fs:6);
                end

                for t = 1:minlen_laser
                    train1 = nonemp_laser{1}{t};
                    train2 = nonemp_laser{2}{t};
                    temp_laser(t,:) = histcounts([train1(train1 >=0 & train1 <= 3);...
                        train2(train2 >=0 & train2 <= 3)+3]',0:1/fs:6);
                end

                results(n).ctrl_perf = exps(ns).ctrl_perf_SPIKE.max_masked{ch}(inds(ii));
                results(n).laser_perf = exps(ns).laser_perf_SPIKE.max_masked{ch}(inds(ii));
                results(n).ctrl_d = exps(ns).ctrl_perf_SPIKE.d_masked{ch}(inds(ii));
                results(n).laser_d = exps(ns).laser_perf_SPIKE.d_masked{ch}(inds(ii));
                results(n).ctrl_driven_FR = calcFRPeriod(spks_ctrl,'driven');
                results(n).laser_driven_FR = calcFRPeriod(spks_laser,'driven');
                results(n).ctrl_spon_FR = calcFRPeriod(spks_ctrl,'laser');
                results(n).laser_spon_FR = calcFRPeriod(spks_laser,'laser');
                results(n).unit = mua;
                results(n).stim = 'Masked';
                results(n).ch = ch;
                results(n).tloc = locs(tloc);
                results(n).mloc = locs(mloc);
                results(n).subject = subjects{ns};
                results(n).num = unitnum;
                results(n).celltype = celltype;
                n = n + 1;
            end
        end
        unitnum = unitnum + 1;
    end
end
unitnum = unitnum - 1;

% implement logic for control and emergent hotspots
d1 = find( [results.ctrl_d] >= d_thresh & [results.ctrl_perf] >= perf_thresh);  % control
d2 = find( ([results.ctrl_d] < d_thresh & [results.ctrl_perf] < perf_thresh) & ...
    ([results.laser_d] >= d_thresh & [results.laser_perf] >= perf_thresh));   % emergent (stringent)

mainset = results(d1);
emergent = results(d2);
singleunits = mainset(strcmp({mainset.unit},'single'));

%% Compare between conditions - Fig S8 for performance and spontaneous FR

fields = {'perf','driven_FR','spon_FR'};
field_titles = {'Performance','Driven FR','Spontaneous FR'};

colors = get(gca,'colororder'); close;
clearvars p_clean p_masked d_clean d_masked

for ff = 1:2 % measures
    clearvars h perf_ctrl_clean perf_laser_clean perf_ctrl_masked perf_laser_masked
    tempstruct = singleunits;
    % figure('units','inches','position',[2 2 1.8*2 2.1*2]);
    figure('units','inches','position',[2 2 2 2.4]);

    for ns = 1:length(subjects)

        clean_inds = strcmp({tempstruct.stim},'Clean') & strcmp({tempstruct.subject},subjects{ns});
        masked_inds = strcmp({tempstruct.stim},'Masked') & strcmp({tempstruct.subject},subjects{ns});

        perf_ctrl_clean{ns} = [tempstruct(clean_inds).(['ctrl_' fields{ff}])];
        perf_laser_clean{ns} = [tempstruct(clean_inds).(['laser_' fields{ff}])];

        perf_ctrl_masked{ns} = [tempstruct(masked_inds).(['ctrl_' fields{ff}])];
        perf_laser_masked{ns} = [tempstruct(masked_inds).(['laser_' fields{ff}])];

        subplot(1,2,1);
        %plot((ones(numel(perf_ctrl_clean{ns}),2).*[1 2])',[perf_ctrl_clean{ns};perf_laser_clean{ns}],'-o','color',colors(ns,:),'markersize',4); hold on;
        plot((ones(numel(perf_ctrl_clean{ns}),2).*[1 2])',[perf_ctrl_clean{ns};perf_laser_clean{ns}],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;
        subplot(1,2,2);
        %h{ns} = plot((ones(numel(perf_ctrl_masked{ns}),2).*[1 2])',[perf_ctrl_masked{ns};perf_laser_masked{ns}],'-o','color',colors(ns,:),'markersize',4); hold on;
        plot((ones(numel(perf_ctrl_masked{ns}),2).*[1 2])',[perf_ctrl_masked{ns};perf_laser_masked{ns}],'-o','color',[0.86 0.86 0.86],'markersize',4); hold on;

    end

    [~,p_clean(ff)] = ttest(horzcat(perf_ctrl_clean{:}),horzcat(perf_laser_clean{:}))
    [~,p_masked(ff)] = ttest(horzcat(perf_ctrl_masked{:}),horzcat(perf_laser_masked{:}))
    [d_clean(ff)] = computeCohen_d(horzcat(perf_ctrl_clean{:}),horzcat(perf_laser_clean{:}),'paired')
    [d_masked(ff)] = computeCohen_d(horzcat(perf_ctrl_masked{:}),horzcat(perf_laser_masked{:}),'paired')
    
    subplot(1,2,1);
    set(gca,'xticklabels',{'Control','Laser'}); sigstar([1 2],p_clean(ff)); set(gca,'fontsize',10);
    set(gca,'xtick',[1 2]); title('Clean','Fontsize',10,'fontweight','normal'); % title(['Clean, p = ' num2str(p_clean)]);
    lim1 = get(gca,'ylim'); xlim([0.5 2.5]);

    if ff == 1, ytickformat('percentage'); end

    subplot(1,2,2);
    set(gca,'xticklabels',{'Control','Laser'});sigstar([1 2],p_masked(ff));  set(gca,'fontsize',10);
    set(gca,'xtick',[1 2],'ytick',[]); title('Masked','Fontsize',10,'fontweight','normal'); %title(['Masked, p = ' num2str(p_masked)]);
    lim2 = get(gca,'ylim'); xlim([0.5 2.5])

    if ff == 1, lims = [50 104];
    else, lims = [lim1 lim2]; end

    % plot means
    subplot(1,2,1); ylim([min(lims) max(lims)]); plot([1 2],mean([horzcat(perf_ctrl_clean{:});horzcat(perf_laser_clean{:})],2,'omitnan'),'-ok','linewidth',2,'markersize',4);
    subplot(1,2,2); ylim([min(lims) max(lims)]); plot([1 2],mean([horzcat(perf_ctrl_masked{:});horzcat(perf_laser_masked{:})],2,'omitnan'),'-ok','linewidth',2,'markersize',4);

    sgtitle([field_titles{ff}],'fontsize',10);

   % h = cellfun(@(x) x(1),h);
    %legend(h,cellfun(@(x) extractBefore(x,'-'),subjects,'uniformoutput',false))

    savefig(gcf,[folder filesep field_titles{ff} ', control SUs'])
    [cc,p] = corrcoef([[tempstruct.ctrl_perf],[tempstruct.laser_perf]],...
       [[tempstruct.(['ctrl_' fields{ff}])],[tempstruct.(['laser_' fields{ff}])]])

end

numSUs = length(unique([singleunits.num]));
numNS = length(unique([singleunits(strcmp({singleunits.celltype},'narrow_spiking')).num]));

save(fullfile(folder,'control_SU_stats.mat'),'p_clean','p_masked','d_clean','d_masked','singleunits');

fprintf('%i SUs with %i NS units: %i clean hotspots, %i masked hotspots\n',numSUs,numNS,length(horzcat(perf_ctrl_clean{:})),length(horzcat(perf_ctrl_masked{:})));

%% Look at responses at hotspots - Figure S7

response_folder = [folder filesep 'Hotspot responses'];
if ~isfolder(response_folder)
    mkdir(response_folder); 
end

locs = [90 45 0 -90];
for n = 1:length(singleunits)
    
    subind = find(strcmp(subjects,singleunits(n).subject));
    ch = singleunits(n).ch;
    tloc = find(locs == singleunits(n).tloc);
    mloc = find(locs == singleunits(n).mloc);

    if isempty(mloc) % clean
        ctrl_spks = exps(subind).ctrl_spks.Spks_clean{ch}(:,tloc,:);
        laser_spks = exps(subind).laser_spks.Spks_clean{ch}(:,tloc,:);
    else
        ctrl_spks = exps(subind).ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:);
        laser_spks = exps(subind).laser_spks.Spks_masked{ch}(:,tloc,mloc,:);
    end

    perfs = [singleunits(n).ctrl_perf singleunits(n).laser_perf];
    lind = [tloc mloc];

    [title_str] = plotCompareResponses(ctrl_spks,laser_spks,0.5,3.5,ch,lind,perfs);
    saveas(gcf,[response_folder filesep subjects{subind} ', ' title_str '.png']);
    savefig(gcf,[response_folder filesep subjects{subind} ', ' title_str]);
    close;

end

%% Plot laser onset PSTHs for NS cells only

% center time vector around laser onset

dt = 2 / 1000;
T_vec = [-fliplr(0.05:dt:0.1) (-0.05+dt):dt:(0+dt)]; 

% pool all NS units and determine if suppression was significant during
% first dt ms of laser

psth_ctrl = []; psth_laser = [];
for ns = 1:length(exps)
    Chans = NS{ns}';
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
plot(T_vec(1:end-1)+0.05,av_ctrl,'k','linewidth',1);

hold on; plot(T_vec(1:end-1)+0.05,av_laser,'r','linewidth',1);
text(0,30,{['p = ' num2str(p)],['d = ' num2str(d)]});

ylabel('FR (Hz)'); xlabel('Time after laser onset (s)'); set(gca,'fontsize',8)
patch([T_vec(1:end-1) fliplr(T_vec(1:end-1))]+0.05,[ll_ctrl fliplr(ul_ctrl)],'black','facealpha',0.3,'edgecolor','none')
patch([T_vec(1:end-1) fliplr(T_vec(1:end-1))]+0.05,[ll_laser fliplr(ul_laser)],'red','facealpha',0.3,'edgecolor','none')

title('Narrow-spiking single-unit activity around laser onset','fontweight','normal');
save([folder filesep 'Control_subject_NS_onset.mat'],'psth_ctrl','psth_laser','T_vec');

savefig(gcf,[folder filesep 'Mean control subject laser onset raster.fig']);

% after: run compareLaserOnsetDiffs (must have also run finalPopAnalysis to get Arch_subject_NS_onset)