% stats for coldspot analysis (see Statistics and Reproducibility section
% in Methods)

load('spikeMetrics.mat')

artifact_chs = {[],29,22,25,[],27,[23 26],[],25};
chanMap = [9 8 6 11 7 10 4 12;
    13 1 15 3 14 2 16 5;
    31 19 29 17 32 20 27 18;
    26 23 21 28 24 25 22 30];

layer_names = {'L1','L2/3','L4','L5','L6'};
layer_edges = [0.5 3.5 6.5 8.5 11.5 15.5];

%% Re-run results but for SUs only

results = struct; n = 1; unitnum = 1;
perf_thresh = 70; d_thresh = 1;

for ns = 1:length(exps)
                
    Chans = SU{ns};

    for ch = Chans
        
        [shank,depth] = find(chanMap == ch);
        temp_dpth = 8 + depth - exps(ns).gran_layers(shank);
        layer_ind = discretize(temp_dpth,layer_edges);
        
        % check if performance is from low # trials
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end
        
        ctrl_res = structfun(@(x) any(x == ch),cellinfo(ns));
        
        if sum(ctrl_res) == 0, celltype = 'onset'; else
        celltype = alltypes{structfun(@(x) any(x == ch),cellinfo(ns))}; end
        
        inds = [];
        
        % clean
        inds = findColdspots(exps(ns).ctrl_perf_SPIKE.Max{ch},exps(ns).ctrl_perf_SPIKE.d{ch},...
            exps(ns).ctrl_spks.Spks_clean{ch},perf_thresh,d_thresh);
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
                results(n).unit = 'SU';
                results(n).stim = 'Clean';
                results(n).ch = ch;
                results(n).tloc = locs(inds(ii));
                results(n).mloc = nan;
                results(n).subject = subjects{ns};
                results(n).num = unitnum;
                results(n).layer = layer_names{layer_ind};
                results(n).type = celltype;
                n = n + 1;
            end
        end
        
        inds = [];
        
        % masked
        inds = findColdspots(exps(ns).ctrl_perf_SPIKE.max_masked{ch},exps(ns).ctrl_perf_SPIKE.d_masked{ch},...
            exps(ns).ctrl_spks.Spks_masked{ch},perf_thresh,d_thresh);
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
                results(n).unit = 'SU';
                results(n).stim = 'Masked';
                results(n).ch = ch;
                results(n).tloc = locs(tloc);
                results(n).mloc = locs(mloc);
                results(n).subject = subjects{ns};
                results(n).num = unitnum;
                results(n).layer = layer_names{layer_ind};
                results(n).type = celltype;
                n = n + 1;
            end
        end
        unitnum = unitnum + 1;
    end
end
unitnum = unitnum - 1;

%% Separate coldspots by effect sizes

d_groups = [0.2 0.5 0.8 30];

% exclude performances below chance value of 50%
better_res = results([results.ctrl_perf] >= 50);

for k = 1:3
    inds = [better_res.ctrl_d] >= d_groups(k) & [better_res.ctrl_d] < d_groups(k+1);
    tempstruct = better_res(inds);
    
    clean_inds = strcmp({tempstruct.stim},'Clean');
    masked_inds = strcmp({tempstruct.stim},'Masked');

    ctrl_clean_perf = [tempstruct(clean_inds).ctrl_perf];
    laser_clean_perf = [tempstruct(clean_inds).laser_perf];
    
    ctrl_masked_perf = [tempstruct(masked_inds).ctrl_perf];
    laser_masked_perf = [tempstruct(masked_inds).laser_perf];

    n_clean(k) = sum(clean_inds);
    n_masked(k) = sum(masked_inds);

    [~,p_clean(k)] = ttest(ctrl_clean_perf,laser_clean_perf);
    [~,p_masked(k)] = ttest(ctrl_masked_perf,laser_masked_perf);

    d_clean(k) = computeCohen_d(ctrl_clean_perf,laser_clean_perf,'paired');
    d_masked(k) = computeCohen_d(ctrl_masked_perf,laser_masked_perf,'paired');

end

% make table
SU_coldspot_stats = table(n_clean',p_clean',d_clean',n_masked',p_masked',d_masked',...
    'VariableNames',{'n_clean','p_clean','d_clean','n_masked','p_masked','d_masked'},...
    'RowNames',{'0.2-0.5','0.5-0.8','>0.8'})

