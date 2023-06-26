function FR = calcFRPeriod(spks,period)

% inputs:
% spks - ntargets x ntrials cell
% period - string ('driven' or 'spon')

[~,nTargets] = size(spks);

% calculate average firing rate across both targets for either spontaneous
% or driven period

for tid = 1:nTargets
    nonemp = find(~cellfun(@isempty,spks(:,tid)));
    if strcmp(period,'driven')
        FR_trials = cellfun(@(x) sum(x >= 0 & x < 3),spks(nonemp,tid))/3;
    elseif strcmp(period,'spontaneous')
        FR_trials = cellfun(@(x) sum(x < -0.05),spks(nonemp,tid))/0.95;
    elseif strcmp(period,'whole')
        FR_trials = cellfun(@numel,spks(nonemp,tid))/5;
    elseif strcmp(period,'laser')
        FR_trials = cellfun(@(x) sum(x >= -0.05 & x < 0),spks(nonemp,tid))/0.05;
    elseif strcmp(period,'prelaser')
        FR_trials = cellfun(@(x) sum(x >= -0.1 & x < -0.05),spks(nonemp,tid))/0.05;
    elseif strcmp(period,'change')
        FR_trials = cellfun(@(x) sum(x >= 0 & x < 3)/3 - sum(x >= -0.05 & x < 0)/0.05,spks(nonemp,tid));
    elseif strcmp(period,'onset')
        FR_trials = cellfun(@(x) sum(x >= 0 & x < 0.5)/0.5,spks(nonemp,tid));
    end
    % mean within each target
    FR(tid) = mean(FR_trials);
end

% mean across all trials
FR = mean(FR);

end