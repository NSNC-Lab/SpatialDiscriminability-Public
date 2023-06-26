function [pt_ratio,pt_int,pt_ch,onset] = calcPeakTrough(exps)
figure;

for ns = 1:length(exps)
    
    % chans_ctrl = findgoodchans(exps(ns).ctrl_perf_SPIKE,perf_thresh,d_thresh);
    % chans_laser = findgoodchans(exps(ns).laser_perf_SPIKE,perf_thresh,d_thresh);
    % Chans = union(chans_ctrl,chans_laser);
    
    Chans = find(~cellfun(@isempty,exps(ns).ctrl_spks.Spks_clean));
    
    if isempty(Chans), continue; end
    
    [pt_ratio{ns},pt_int{ns},pt_ch{ns}] = calcPTMetrics(exps(ns).ctrl_spks.snips_clean,exps(ns).ctrl_spks.snips_masked,Chans);
    [onset{ns}] = calcOnsetResponse(exps(ns).ctrl_spks.Spks_masked,Chans);
    
    % scatter(pt_int{ns},onset{ns},'filled','k'); hold on;
end

% xlabel('PT interval (ms)'); ylabel('Onset response (ms)');
% xlim([0.2 1]); %ylim([0 50]);

end

function [pt_ratio,pt_int,pt_ch] = calcPTMetrics(snips_clean,snips_masked,Chans)

n = 1;

for ch = Chans
    
    snips = vertcat(snips_clean{ch}{:},snips_masked{ch}{:});
    
    pt_ch(n) = ch;
    
    winlen = size(snips,2);
    spktimes = winlen/4;
    
    %     [peak,peak_ind] = max(snips(:,spktimes:end) - mean(snips,2),[],2);
    %     [trough,trough_ind] = min(snips(:,spktimes-2:spktimes+2) - mean(snips,2),[],2);
    %     pt_int(n) = mean(1000*(peak_ind - trough_ind - 2) / 24414.0625);
    %     pt_ratio(n) = mean(peak./abs(trough));
    
    mean_snip = mean(snips);
    [peak,peak_ind] = max(mean_snip(spktimes:end));
    [trough,trough_ind] = min(mean_snip(spktimes-2:spktimes+2));
    pt_int(n) = 1000*(peak_ind - trough_ind - 2) / 24414.0625;
    pt_ratio(n) = peak./abs(trough);
    
    
    
    
    % scatter(pt_int(ch),pt_ratio(ch),'filled','k');
    
    n = n+1;
end

end


function onset = calcOnsetResponse(spks,Chans)

T_vec = -0.05:1/250:0.05;

binsize = mean(diff(T_vec));
T_vec(end+1) = T_vec(end) + binsize;

n = 1;
for ch = Chans
        
    psth = [];
    for t = 1:numel(spks{ch})
        psth = cat(1,psth,histcounts(spks{ch}{t},T_vec));
    end
    
    % look only at 
    
    psth = sum(psth);
    
    [~,ind] = max(psth);
    onset(n) = max(0,T_vec(ind)*1000);
    
    n = n+1;
end

end