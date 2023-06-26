function title_str = plotCompareResponses(spks1,spks2,cond1,cond2,T_before,T_after,ch,locs,perfs)

% Inputs:
% [ctrl or laser]_spks: Spks_clean or Spks_masked
% T_before and T_after: limits of x-axis
% ch: channel
% locs: 2-element vector of configuration info
% perfs: 2-element vector of control and laser performances
% folder: where to save

spks1 = squeeze(spks1); spks2 = squeeze(spks2);

perfs = round(perfs);

if numel(locs) == 1, locs(2) = NaN; end
[tdeg,mdeg] = getDegrees(locs);

[soundstim(:,1),Fs] = audioread('200k_target1.wav');
soundstim(:,2) = audioread('200k_target2.wav');
title_str = ['CH' num2str(ch) ', T' tdeg char(176)];

% load masked stimulus for masked trials
if ~isnan(locs(2)), masker = audioread('200k_masker1.wav');
    soundstim = soundstim + masker;
    title_str = ['CH' num2str(ch) ', T' tdeg char(176) ', M' mdeg char(176)];
end

Fs = Fs/4; soundstim = downsample(soundstim,4);

nzb = sum(-T_before:1/Fs:T_after < 0);
nza = sum(-T_before:1/Fs:T_after > size(soundstim,1)/Fs);

soundstim = [zeros(nzb,2) ; soundstim ; zeros(nza,2)];

t_stim = linspace(-T_before,T_after,size(soundstim,1));

bin = 0.02; % in seconds

figure('unit','inches','position',[7 6 5.5 4])

% annotate with config
annotation('textbox',[0.58 0.82 0.3 0.1],'String',title_str,'edgecolor','none','fontsize',10)
calcAndAnnotFR(spks1,cond1,[0.58 0.77 0.3 0.1])
calcAndAnnotFR(spks2,cond2,[0.58 0.72 0.3 0.1])


% performance barplot
subplot(5,2,1); hold on;
origpos = get(gca,'position');
origpos(1) = origpos(1) + origpos(3)/2;
origpos(3) = origpos(3)/2; set(gca,'position',origpos)

bar(1,perfs(1),'facecolor','none','linewidth',2,'EdgeColor','k');
yperf = -10;
if perfs(1) < 55, yperf = 10; end
text(0.75,perfs(1)+yperf,[num2str(perfs(1)) '%'])

bar(2,perfs(2),'facecolor','none','linewidth',2,'EdgeColor','r');
if perfs(2) < 55, yperf = 10; end
text(1.75,perfs(2)+yperf,[num2str(perfs(2)) '%'])

axis tight; box off; xlim([0.6 2.4])
ylim([40 max(perfs)+5]);
title('Performance','fontweight','normal'); ytickformat('percentage'); set(gca,'xticklabels',{cond1,cond2},'xtick',[1 2])

for tid = 1:2 % per target

    % stimulus
    subplot(5,2,tid+2);
    plot_stimulus(t_stim,soundstim(:,tid),T_before,T_after)

    % control raster
    subplot(5,2,tid+4);
    nonemp = ~cellfun(@isempty,spks1(:,tid));
    plotSpikeRasterFs(spks1(nonemp,tid)','PlotType','vertline');
    axis tight; xlim([-T_before T_after]); box off; set(gca,'xtick',[])
    ylabel(cond1);

    % laser raster
    subplot(5,2,tid+6);
    LineFormat = struct(); LineFormat.Color = [1 0 0];
    nonemp = ~cellfun(@isempty,spks2(:,tid));
    plotSpikeRasterFs(spks2(nonemp,tid)','PlotType','vertline','LineFormat',LineFormat);
    
    axis tight; xlim([-T_before T_after]); box off; set(gca,'xtick',[])
    ylabel(cond2);

    % control and laser psth
    subplot(5,2,tid+8)
    makePSTH(spks1(:,tid),-T_before:bin:T_after,bin,'k'); hold on;
    makePSTH(spks2(:,tid),-T_before:bin:T_after,bin,'r')

    ylims(tid,:) = get(gca,'ylim');
    axis tight; xlim([-T_before T_after]); box off
    ylabel('FR [Hz]'); xlabel('Time [s]');

end

% adjust PSTH y-lims

ymax = max(ylims,[],'all');
subplot(5,2,9); ylim([0 ymax]);
subplot(5,2,10); ylim([0 ymax]);

end

function plot_stimulus(tid,y,T_before,T_after)

plot(tid,y,'k')
set(gca,'XTick', [],'Ytick',[]);
ylim([-1 1])
xlim([-T_before T_after]);
axis tight; box off;

end

function [tdegree,mdegree] = getDegrees(locs)

switch locs(1)
    case 1, tdegree='90';
    case 2, tdegree='45';
    case 3, tdegree='0';
    case 4, tdegree='-90';
end

switch locs(2)
    case 1, mdegree='90';
    case 2, mdegree='45';
    case 3, mdegree='0';
    case 4, mdegree='-90';
    otherwise, mdegree=NaN;
end

end

function [ll,ul,av,psth] = makePSTH(spks,T_vec,bin,color)

psth = [];
for n = 1:numel(spks)
    psth = cat(1,psth,histcounts(spks{n},T_vec));
end

[ll,ul,av] = sem(psth);
plot(T_vec(1:end-1),av/bin,color,'linewidth',1);
patch([T_vec(1:end-1) fliplr(T_vec(1:end-1))],[ll fliplr(ul)]/bin,color,'facealpha',0.3,'edgecolor','none');

end

function calcAndAnnotFR(spks,cond,position)

temp = cellfun(@(x) numel(x >= 0 & x < 3),spks(:));
temp(temp == 0) = [];
temp = temp/3;
fr = mean(temp); fr_std = std(temp);

frStr = sprintf('%s FR: %0.2f ± %0.2f Hz',cond,fr,fr_std);

annotation('textbox',position,...
    'string',frStr,...
    'FitBoxToText','on',...
    'LineStyle','none',...
    'fontsize',10)

end

