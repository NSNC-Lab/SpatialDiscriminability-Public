function [TS,RMS] = calcTrialSim(spks)

T_vec = 0:0.025:3;

% trial similarity
for t = 1:2
    
    temp = [];
    for n = 1:100 %resample trials 100 times
        inds = randperm(numel(spks(:,t)));
        inds1 = inds(1:floor(numel(inds)/2)); inds2 = setdiff(inds,inds1);
        
        PSTH1 = histcounts(vertcat(spks{inds1,t}),T_vec);
        PSTH2 = histcounts(vertcat(spks{inds2,t}),T_vec);
        
        c = corrcoef(PSTH1,PSTH2);
        
        temp = cat(1,temp,c(1,2));
        
        % trial-average PSTHs for RMS difference
    end
    TS(t) = mean(temp);
end

TS = mean(TS);

for t = 1:2
    temp = sum(~cellfun(@isempty,spks(:,t)));
    PSTH(t,:) = histcounts(vertcat(spks{:,t}),T_vec) / temp; %/ 0.025; % trial-average and convert to Hz
    
    %normalize PSTH so that sum is 1
    PSTH(t,:) = PSTH(t,:) / sum(PSTH(t,:));
end

RMS = rms(diff(PSTH));

end