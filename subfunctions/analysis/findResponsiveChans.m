function Chans = findResponsiveChans(ctrl_spks)

allChans = find(~cellfun(@isempty,ctrl_spks.Spks_clean));
Chans = [];

for ch = allChans
    
    % average firing rate
    
    clean_spks = ctrl_spks.Spks_clean{ch};
    masked_spks = ctrl_spks.Spks_masked{ch};
    FR_0 = []; p1 = []; p2 = [];
    for n = 1:4
        spks = squeeze(clean_spks(:,n,:));
        % FR_0 = cat(1,FR_0,cellfun(@(x) sum(x < 0),spks(:)));
        FR_0 = cellfun(@(x) sum(x < 0),spks(:));
        p1 = cat(2,p1,maxFRResp(spks,FR_0));
        p2 = cat(2,p2,avgFRResp(spks,FR_0));
    end
    
    for n = 1:16
        [tloc,mloc] = ind2sub([4 4],n);
        spks = squeeze(masked_spks(:,tloc,mloc,:));
        % FR_0 = cat(1,FR_0,cellfun(@(x) sum(x < 0),spks(:)));
        FR_0 = cellfun(@(x) sum(x < 0),spks(:));
        p1 = cat(2,p1,maxFRResp(spks,FR_0));
        p2 = cat(2,p2,avgFRResp(spks,FR_0));
    end        
    
    % if any(p1 < 0.05) || any(p2 < 0.05)
    if sum(p1 < 0.05)>20 || sum(p2 < 0.05)>20
        Chans = cat(2,Chans,ch);
    end
    
end

end

function p = maxFRResp(spks,FR_0)

T_vec = -1:1/50:4;

for t = 1:2
    
    spkset = spks(:,t);
    nonemp = sum(~cellfun(@isempty,spkset));
    temp = cellfun(@(x) histcounts(x,T_vec),spkset,'uniformoutput',false);
    
    PSTH = sum(vertcat(temp{:}))  / nonemp / (1/50);
    
    % find maximal firing rate
    inds = find(T_vec(1:end-1) >= 0 & T_vec(1:end-1) < 3);
    [~,maxind] = max(PSTH(inds));
    best_time = T_vec(inds(maxind));
    
    % estimate FR
    max_resp = cellfun(@(x) sum(x >= best_time & x < best_time + 1/50),spkset) / (1/50);
    
    [~,p(t)] = ttest2(FR_0,max_resp','vartype','unequal','alpha',0.05);
    
end

end

function p = avgFRResp(spks,FR_0)

for t = 1:2
    
    spkset = spks(:,t);
    temp = cellfun(@(x) sum(x >= 0 & x < 3)/3,spkset(~cellfun(@isempty,spkset)));
        
    % estimate FR
    [~,p(t)] = ttest2(FR_0,temp,'vartype','unequal','alpha',0.05);
    
end

end