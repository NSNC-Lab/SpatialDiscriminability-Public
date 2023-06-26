function [results,perf_comps] = popSearchHotspots(exps,subjects,perf_thresh,d_thresh,cellinfo,SU,artifact_chs)

locs = [90 45 0 -90];
layer_names = {'L1','L2/3','L4','L5','L6'};
layer_edges = [0.5 3.5 6.5 8.5 11.5 15.5];

chanMap = [9 8 6 11 7 10 4 12;
    13 1 15 3 14 2 16 5;
    31 19 29 17 32 20 27 18;
    26 23 21 28 24 25 22 30];

results = struct; n = 1; unitnum = 1;

alltypes = fieldnames(cellinfo);

fs = 24414.0625;

for ns = 1:length(exps)
                
    Chans = findResponsiveChans(exps(ns).ctrl_spks);
    Chans = setdiff(Chans,artifact_chs{ns});
    
    for ch = Chans
        
        [shank,depth] = find(chanMap == ch);
        temp_dpth = 8 + depth - exps(ns).gran_layers(shank);
        layer_ind = discretize(temp_dpth,layer_edges);
                
        if ismember(ch,SU{ns})
            mua = 'single';
        else
            mua = 'mua';
        end
        
        % check if performance is from low # trials
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end
        
        ctrl_res = structfun(@(x) any(x == ch),cellinfo(ns));
        
        if sum(ctrl_res) == 0, celltype = 'reg_spiking'; else
        celltype = alltypes{structfun(@(x) any(x == ch),cellinfo(ns))}; end
        
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
                results(n).ctrl_tau = fineSearchTau(spks_ctrl,exps(ns).ctrl_perf.opt_tau{ch}(inds(ii)));
                results(n).laser_tau = fineSearchTau(spks_laser,exps(ns).laser_perf.opt_tau{ch}(inds(ii)));
                [results(n).ctrl_TS,results(n).ctrl_RMS] = calcTrialSim(spks_ctrl);
                [results(n).laser_TS,results(n).laser_RMS] = calcTrialSim(spks_laser);
                results(n).unit = mua;
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
                results(n).ctrl_tau = fineSearchTau(spks_ctrl,exps(ns).ctrl_perf.opt_tau_masked{ch}(inds(ii)));
                results(n).laser_tau = fineSearchTau(spks_laser, exps(ns).laser_perf.opt_tau_masked{ch}(inds(ii)));
                [results(n).ctrl_TS,results(n).ctrl_RMS] = calcTrialSim(spks_ctrl);
                [results(n).laser_TS,results(n).laser_RMS] = calcTrialSim(spks_laser);
                results(n).unit = mua;
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

%% For different performance metrics
% 
perf_comps = struct;

n = 1;

% distance matrix example nspots = 1 (clean)
for ns = 1:length(exps)
    
    Chans = findResponsiveChans(exps(ns).ctrl_spks);
    Chans = setdiff(Chans,artifact_chs{ns});

    for ch = Chans
        
        % check if performance is from low # trials
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end

        [shank,depth] = find(chanMap == ch);
        temp_dpth = depth - exps(ns).gran_layers(shank) + 8;
        layer = discretize(temp_dpth,layer_edges);
        
        if ismember(ch,SU{ns})
            mua = 'single';
        else
            mua = 'mua';
        end
               
        temp1 = structfun(@(x) any(x == ch),cellinfo(ns));
        if sum(temp1) == 0, celltype = 'reg_spiking'; else
            celltype = alltypes{structfun(@(x) any(x == ch),cellinfo(ns))}; end
        
        % clean
        inds_ctrl = findHotspots(exps(ns).ctrl_perf_SPIKE.Max{ch},exps(ns).ctrl_perf_SPIKE.d{ch},...
            exps(ns).ctrl_spks.Spks_clean{ch},perf_thresh,d_thresh);
        inds_laser = findHotspots(exps(ns).laser_perf_SPIKE.Max{ch},exps(ns).laser_perf_SPIKE.d{ch},...
            exps(ns).laser_spks.Spks_clean{ch},perf_thresh,d_thresh);
        inds = union(inds_ctrl,inds_laser);
        if isrow(inds), inds = inds'; end

        if isrow(inds)
            inds = inds';
        end
        
        if ~isempty(inds)
            for ii = 1:length(inds)
                
                perf_comps(n).ctrl_perf_SPIKE = exps(ns).ctrl_perf_SPIKE.Max{ch}(inds(ii));
                perf_comps(n).laser_perf_SPIKE = exps(ns).laser_perf_SPIKE.Max{ch}(inds(ii));
                perf_comps(n).ctrl_perf_RISPIKE = exps(ns).ctrl_perf_RI_SPIKE.Max{ch}(inds(ii));
                perf_comps(n).laser_perf_RISPIKE = exps(ns).laser_perf_RI_SPIKE.Max{ch}(inds(ii));
                perf_comps(n).ctrl_perf_ISI = exps(ns).ctrl_perf_ISI.Max{ch}(inds(ii));
                perf_comps(n).laser_perf_ISI = exps(ns).laser_perf_ISI.Max{ch}(inds(ii));
                perf_comps(n).ctrl_perf_spkct = exps(ns).ctrl_perf_spkct.Max{ch}(inds(ii));
                perf_comps(n).laser_perf_spkct = exps(ns).laser_perf_spkct.Max{ch}(inds(ii));
                perf_comps(n).ctrl_d = exps(ns).ctrl_perf_SPIKE.d{ch}(inds(ii));
                perf_comps(n).laser_d = exps(ns).ctrl_perf_SPIKE.d{ch}(inds(ii));
                perf_comps(n).unit = mua;
                perf_comps(n).ch = ch;
                perf_comps(n).subject = subjects{ns};
                perf_comps(n).stim = 'Clean';
                perf_comps(n).tloc = locs(inds(ii));
                perf_comps(n).mloc = nan;
                perf_comps(n).layer = layer;
                perf_comps(n).type = celltype;
                n = n + 1;
            end
        end
        
        % masked
        inds_ctrl = findHotspots(exps(ns).ctrl_perf_SPIKE.max_masked{ch},exps(ns).ctrl_perf_SPIKE.d_masked{ch},...
            exps(ns).ctrl_spks.Spks_masked{ch},perf_thresh,d_thresh);
        inds_laser = findHotspots(exps(ns).laser_perf_SPIKE.max_masked{ch},exps(ns).laser_perf_SPIKE.d_masked{ch},...
            exps(ns).laser_spks.Spks_masked{ch},perf_thresh,d_thresh);
        inds = union(inds_ctrl,inds_laser);
        if isrow(inds), inds = inds'; end
        
        if isrow(inds)
            inds = inds';
        end
        
        if ~isempty(inds)
            for ii = 1:length(inds)
                [mloc,tloc] = ind2sub([4 4],inds(ii));
                mloc_spk = 5 - mloc;
                
                perf_comps(n).ctrl_perf_SPIKE = exps(ns).ctrl_perf_SPIKE.max_masked{ch}(inds(ii));
                perf_comps(n).laser_perf_SPIKE = exps(ns).laser_perf_SPIKE.max_masked{ch}(inds(ii));
                perf_comps(n).ctrl_perf_RISPIKE = exps(ns).ctrl_perf_RI_SPIKE.max_masked{ch}(inds(ii));
                perf_comps(n).laser_perf_RISPIKE = exps(ns).laser_perf_RI_SPIKE.max_masked{ch}(inds(ii));
                perf_comps(n).ctrl_perf_ISI = exps(ns).ctrl_perf_ISI.max_masked{ch}(inds(ii));
                perf_comps(n).laser_perf_ISI = exps(ns).laser_perf_ISI.max_masked{ch}(inds(ii));
                perf_comps(n).ctrl_perf_spkct = exps(ns).ctrl_perf_spkct.max_masked{ch}(inds(ii));
                perf_comps(n).laser_perf_spkct = exps(ns).laser_perf_spkct.max_masked{ch}(inds(ii));
                perf_comps(n).ctrl_d = exps(ns).ctrl_perf_SPIKE.d_masked{ch}(inds(ii));
                perf_comps(n).laser_d = exps(ns).ctrl_perf_SPIKE.d_masked{ch}(inds(ii));
                perf_comps(n).unit = mua;
                perf_comps(n).ch = ch;
                perf_comps(n).subject = subjects{ns};
                perf_comps(n).stim = 'Masked';
                perf_comps(n).tloc = locs(tloc);
                perf_comps(n).mloc = locs(mloc_spk);
                perf_comps(n).layer = layer;
                perf_comps(n).type = celltype;
                n = n + 1;
            end
        end
    end
end

end