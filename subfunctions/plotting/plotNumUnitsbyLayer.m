function plotNumUnitsbyLayer(exps,SU,folder,artifact_chs)


figure('unit','inches','position',[3 3 2.5 2]); hold on;
layer_edges = [0.5 3.5 6.5 8.5 11.5 15.5];
chanMap = [9 8 6 11 7 10 4 12;
    13 1 15 3 14 2 16 5;
    31 19 29 17 32 20 27 18;
    26 23 21 28 24 25 22 30];

SU_counts = zeros(1,5);
MU_counts = zeros(1,5);

for ns = 1:9
    
    Chans = findResponsiveChans(exps(ns).ctrl_spks);
    Chans = setdiff(Chans,artifact_chs{ns});
    
    for ch = Chans
        
        % exclude channels with low trial count
        if any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).ctrl_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        elseif any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_clean{ch}) >= 3,'all') || ...
                any(cellfun(@(x) sum(x == 0),exps(ns).laser_spks.num_spks_masked{ch}) >= 3,'all')
            continue
        end

        [shank,depth] = find(chanMap == ch);
        temp_dpth = 8 + depth - exps(ns).gran_layers(shank);
        layer_ind = discretize(temp_dpth,layer_edges);
                
        if ~ismember(ch,SU{ns})
            MU_counts(layer_ind) = MU_counts(layer_ind) + 1;
        else % singleunit
            SU_counts(layer_ind) = SU_counts(layer_ind) + 1;
        end
    end
    
end

bar([1:5]-0.2,SU_counts,0.4,'facecolor','none','linewidth',2); hold on;
bar([1:5]+0.2,MU_counts,0.4,'facecolor','k','linewidth',2); 
ylabel('# units'); xlabel('Layer');
set(gca,'xticklabels',{'L1','L2/3','L4','L5','L6'},'xtick',1:5,'fontsize',8);
xlim([0.3 5.7]);
legend('SU','MU','location','best');

saveas(gcf,[ folder filesep '# units by layer.png']);
savefig(gcf,[ folder filesep '# units by layer']);
% close;

end