function opt = fineSearchTau(spks,tau)

all_taus = 0.001 * [1,2,4,8,16,32,64,128,256];

orig = find(all_taus == tau);
if orig == 1
    opt = all_taus(orig); return
elseif orig == length(all_taus)
    Tau = all_taus(end-1) : 1/1000: all_taus(end);
else
    Tau = all_taus(orig - 1) : 1/1000: all_taus(orig+1);
end

for t = 1:2
    nonemp{t}= find(~cellfun(@isempty,spks(:,t)));
    len(t) = length(nonemp{t});
end

for t = 1:2 %stimuli#
    for k = 1:min(len)
        temp = spks{nonemp{t}(k)};
        spTimes{k,t} = temp(temp >= 0 & temp < 3);
    end
end

distMat = calcvr(spTimes,Tau);
perf = calcpcStatic(distMat, k, 2, 0);

[~,multi] = max(perf);

opt = Tau(multi);
    
end