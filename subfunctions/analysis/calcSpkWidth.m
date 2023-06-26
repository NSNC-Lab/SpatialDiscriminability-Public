function [spk_width] = calcSpkWidth(snips_clean,snips_masked)

snips = vertcat(snips_clean{:},snips_masked{:});

winlen = size(snips,2);
spktimes = winlen/4; %1/4 time window of spike snippet is usually where spike time is assigned

mean_snip = mean(snips);
[~,peak_ind] = max(mean_snip(spktimes-2:end));
[~,trough_ind] = min(mean_snip(spktimes-2:end));
spk_width = 1000*(peak_ind - trough_ind) / 24414.0625;

end
