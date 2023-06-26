function maxFR = estimateMaxFR(spks)

T_vec = -1:1/50:4;
for t = 1:2

    spkset = spks(:,t);
    nonemp = sum(~cellfun(@isempty,spkset));
    temp = cellfun(@(x) histcounts(x,T_vec),spkset,'uniformoutput',false);
    
    PSTH = sum(vertcat(temp{:})) * (50 / nonemp);
    
    inds = find(T_vec(1:end-1) >= 0 & T_vec(1:end-1) < 3);
    maxFR(t) = max(PSTH(inds));
    
end

maxFR = max(maxFR);

end