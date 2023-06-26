function [mean_onset] = plotOnsetDepth(exps,power,TMR)

shnks = fieldnames(exps(1).csd_onset);
tempdata = [];

for ns = 1:length(exps)
                
    onsets = zeros(15,4);
    
    % granular layer will be in the middle, i.e. row 8
    for g = 1:4
        inds = (2:7) + (8 - exps(ns).gran_layers(g));
        onsets(inds,g) = exps(ns).csd_onset.(shnks{g});
    end
    
    tempdata = cat(2,tempdata,onsets);
    
end

tempdata = tempdata';

figure('unit','inches','position',[4 4 2.25 2]);

depth = -700:100:700;

% don't include dummy data
onset_lim = 0.5;
[rows,depth_inds] = find(tempdata ~= 0 & tempdata < onset_lim);

for i = 1:length(rows)
     good_data(i) = tempdata(rows(i),depth_inds(i));
end

scatter(depth(depth_inds),1000*good_data,50,'k','linewidth',1); hold on;
for d = 1:length(depth)
    mean_onset(d) = mean(good_data(depth_inds == d));
end
plot(depth,1000*mean_onset,'k','linewidth',2);
xlim([-600 600]); xlabel('Depth ({\mu}m)');
ylim([0 onset_lim*1000]);
ylabel('CSD sink onset (ms)');

title('CSD sink onset vs. depth');
saveas(gcf,['CSD sink onset vs. depth, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR.png']);
savefig(gcf,['CSD sink onset vs. depth, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR.fig']);

xl = table(depth(depth_inds)',good_data','VariableNames',{'DepthIndex','Onset'})

writetable(xl,['CSD sink onset vs. depth, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR.xls']);

end
