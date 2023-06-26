function plot_all_grids(TMR,power,subfolder,type,d1,duds,cluster_chs,cluster_ids)

%%% INPUTS
% TMR: target-to-masker ratio
% figuresdir: subject name
% type: 'control' or 'laser'
% d1: struct of spiking data, should be the same as type

% check if faulty channels are removed
chFlag = length(d1.max_masked) == (32 - length(duds));

nskips = 0;

for ch = 1:length(d1.max_masked)
    
    if ~isempty(duds)
        if chFlag == 0 && ismember(ch,duds) % if dud channel isn't taken out
            continue
        elseif chFlag == 1 && ch >= duds(1)   % if dud channel is taken out
            nskips = sum(ch >= duds);
        end
    end
    
    if isempty(d1.Max{ch})
        continue
    end
        
    grid_folder = [subfolder,'/Grids'];
    
    if ~isfolder(grid_folder)
        mkdir(grid_folder);
    end
    
    barmin = 50;
    barmax = 90;
    
    plot_spatial_grid(ch,ch+nskips,type,TMR,d1.Max{1,ch},d1.max_masked{1,ch},barmin,barmax)
        
    if exist('cluster_chs','var')
        
        h = findobj(gcf,'type','axes');  
        axes(h(2));
        title(sprintf('%s Cluster %g (CH%g), %gdB',type,cluster_ids(ch),cluster_chs(ch),TMR));

        saveas(gcf,[grid_folder '/' type num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE_V2.png']);
        savefig(gcf,[grid_folder '/' type num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE_V2.fig']);
    else
        saveas(gcf,[grid_folder '/' type ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE_V2.png']); 
        savefig(gcf,[grid_folder '/' type ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE_V2.fig']);
    end
    
    close;
end

end