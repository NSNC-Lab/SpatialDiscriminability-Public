function plotSpikeWidthvsFR(exps,SU)

spon_FR = []; spk_width = []; subs = []; % chs = [];
for ns = 1:length(exps)
    
    Chans = SU{ns};

    for ch = Chans
        
        % exclude channels with low trial count
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end

        spk_width = cat(1,spk_width,calcSpkWidth(exps(ns).ctrl_spks.snips_clean{ch},exps(ns).ctrl_spks.snips_masked{ch}));
        
        temp = [];
        for t = 1:4
            temp = cat(1,temp,squeeze(calcFRPeriod(exps(ns).ctrl_spks.Spks_clean{ch}(:,t,:),'spontaneous')));
        end
        for t = 1:16
            [tloc,mloc] = ind2sub([4 4],t);
            temp = cat(1,temp,squeeze(calcFRPeriod(exps(ns).ctrl_spks.Spks_masked{ch}(:,tloc,mloc,:),'spontaneous')));
        end
        spon_FR = cat(1,spon_FR,mean(temp));
        subs = cat(1,subs,ns*ones(numel(mean(temp)),1));
        % chs = cat(1,chs,Chans);
    end
    
end

ns_inds = spk_width < 0.5;

figure('unit','inches','position',[3 3 3 2.5]);

scatter(spk_width(~ns_inds),spon_FR(~ns_inds),25,[0 0.4470 0.7410],'^','k'); hold on; % regular spiking
scatter(spk_width(ns_inds),spon_FR(ns_inds),25,[0.8500 0.3250 0.0980],'o','k');  % narrow-spiking

xlabel('Spike width (ms)'); ylabel('Spontaneous FR (Hz)'); set(gca,'fontsize',8);

% plot mean scatters

scatter(mean(spk_width(~ns_inds)),mean(spon_FR(~ns_inds)),75,'^','k','filled');
scatter(mean(spk_width(ns_inds)),mean(spon_FR(ns_inds)),75,'o','k','filled');

[yll,yul,yav] = sem(spon_FR(~ns_inds));
[xll,xul,xav] = sem(spk_width(~ns_inds));

errorbar(mean(spk_width(~ns_inds)),mean(spon_FR(~ns_inds)),yav-yll,yul-yav,[],[],'k','linewidth',1);

[yll,yul,yav] = sem(spon_FR(ns_inds));
[xll,xul,xav] = sem(spk_width(ns_inds));

errorbar(mean(spk_width(ns_inds)),mean(spon_FR(ns_inds)),yav-yll,yul-yav,[],[],'k','linewidth',1);

[~,p] = ttest2(spon_FR(~ns_inds),spon_FR(ns_inds),'vartype','unequal')

[d] = computeCohen_d(spon_FR(~ns_inds),spon_FR(ns_inds),'independent')

line([0.5 0.5],[0 40],'linestyle','--','color','k')
legend('RS','NS','location','best'); xlim([0 1]);

%saveas(gcf,'Population results/Spike width vs spontaneous FR.png');
%savefig(gcf,'Population results/Spike width vs spontaneous FR');

end

function [spk_width] = calcSpkWidth(snips_clean,snips_masked)

snips = vertcat(snips_clean{:},snips_masked{:});

winlen = size(snips,2);
spktimes = winlen/4; %1/4 time window of spike snippet is usually where spike time is assigned

mean_snip = mean(snips);
[~,peak_ind] = max(mean_snip(spktimes-2:end));
[~,trough_ind] = min(mean_snip(spktimes-2:end));
spk_width = 1000*(peak_ind - trough_ind) / 24414.0625;

end
