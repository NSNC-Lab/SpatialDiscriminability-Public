function trialinfo_new = restoreTS(trialinfo)

% col 1: target location
% col 2: masker location
% col 3: target identity
% col 4: laser on/off

% we'll do a running counter for the # trials per type later

trialinfo_new = trialinfo;

trialinfo_new(trialinfo_new(:,2) == 16,2) = 0;

for n = 1:size(trialinfo,1)
    if trialinfo_new(n,3) > 10 && trialinfo_new(n,3) <= 20
        trialinfo_new(n,3) = 1;
        trialinfo_new(n,2) = trialinfo_new(n,1);
    elseif trialinfo_new(n,3) > 20 && trialinfo_new(n,3) <= 30
        trialinfo_new(n,3) = 2;
        trialinfo_new(n,2) = trialinfo_new(n,1);
    end
end

trialinfo_new(:,4) = trialinfo_new(:,5);

end