function plotRastersTarget(in1,N_cum1,in2,N_cum2,type1,chann,T_before,T_after,...
    T_vec,Max,subfolder)
% in1 is spike times, i.e Spks_clean{1,ch}, loc is which location i.e [speaker,song]
% from {x,z,k} = {speaker,song,trial} so c = x
% N_cum1 is the binned spike count n_cum_clean{1,ch}{speaker,song}
% (cum is summed across all trials)
% type1 is title (laser/control)
% chan is channel Max is max number of spikes in the binned
% in1,sti1,n_cum1 == MASKERS, the 2's are the clean

% plotRastersTarget(Spks_masked,n_cum_masked,Spks_clean,...
%    n_cum_clean,'control',[1:32],1,4,-1:1/50:4,50,sub)

if ~isfolder([subfolder,'/Rasters'])
    mkdir([subfolder,'/Rasters'])
end

% call figure here for efficiency
width = 8.5; height = 4.5;
figure('visible','off','unit','inches','position',[2 2 width height]);
% figuresize(width,height,gcf,'Inches');
% set(gcf,'PaperUnits','inches');

if isempty(chann)
   return 
end

% load sound stimuli first
[soundstim(:,1),Fs] = audioread('200k_target1.wav');
soundstim(:,2) = audioread('200k_target2.wav');
soundstim(:,3) = audioread('200k_masker1.wav');

Fs = Fs/4;
soundstim = downsample(soundstim,4);

stimvec = -T_before:1/Fs:T_after;

nzb = sum(stimvec < 0);
nza = sum(stimvec > size(soundstim,1)/Fs);
            
soundstim = [zeros(nzb,3);soundstim;zeros(nza,3)];
                    
t_stim = linspace(-T_before,T_after,size(soundstim,1));

% define subplot heights and width
w_stim = 0.42;
h_stim = 0.166;
h_rast = 0.25;

h_space = 0.04;
w_space = 0.04;

y0 = 0.15;
x0 = 0.08;

for ch = 1:length(chann)
            
    disp([type1,' CH',num2str(chann(ch))]);
    
    if ~isfolder([subfolder,'/Rasters/CH',num2str(chann(ch))])
        mkdir([subfolder,'/Rasters/CH',num2str(chann(ch))])
    else
        continue
    end
    
    % MASKED RASTERS
    
    loc = [];
    for i = 1:4
        for j = 1:4
            for z = 1:2
                
                loc(1) = i;   % target location
                loc(2) = j;   % masker location
                loc(3) = z;   % target identity
                
                switch i
                    case 1
                        tdegree='90';
                    case 2
                        tdegree='45';
                    case 3
                        tdegree='0';
                    case 4
                        tdegree='-90';
                end
                
                switch j
                    case 1
                        mdegree='90';
                    case 2
                        mdegree='45';
                    case 3
                        mdegree='0';
                    case 4
                        mdegree='-90';
                end
                                
                spikes{5} = [];
                nonemp = find(~cellfun(@isempty,in1{1,chann(ch)}(:,i,j,z)));
                
                if length(nonemp) <= 3 % for TM configurations with no data
                    break
                end
                
                for x=1:length(nonemp)
                    if size(in1{1,chann(ch)}{nonemp(x),loc(1),loc(2),loc(3)},1) > size(in1{1,chann(ch)}{nonemp(x),loc(1),loc(2),loc(3)},2)
                        spikes{x} = in1{1,chann(ch)}{nonemp(x),loc(1),loc(2),loc(3)}';
                    else
                        spikes{x} = in1{1,chann(ch)}{nonemp(x),loc(1),loc(2),loc(3)};
                    end
                end
                            
                % plot stimulus on top
                subplot('position',...
                    [x0+(w_stim+w_space)*(z-1) y0+(h_rast*2 +h_space*2) w_stim h_stim]);
                axis tight
                
                plot(t_stim,soundstim(:,loc(3))+soundstim(:,3),'k')
                set(gca,'XTick', [],'ytick',[]);
                ylim([-1 1])
                xlim([-T_before T_after]); box off
                
                hold on
                
                % plot raster in middle
                subplot('position',...
                    [x0+(w_stim+w_space)*(z-1) y0+(h_rast+h_space) w_stim h_rast]);
                
                for sp = 1:length(nonemp)
                    plot([spikes{sp};spikes{sp}],[ones(size(spikes{sp}))*sp;ones(size(spikes{sp}))*sp-1],'-k','linewidth',.2)
                    hold on
                end
                
                ylim([0 length(nonemp)])
                xlim([-T_before T_after])
                set(gca,'TickDir','out') % draw the tick marks on the outside
                set(gca,'XTick', []) % don't draw y-axis ticks
                %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
                %set(gca,'Color',get(gcf,'Color')) % match figure background
                %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
                box off
                set(gca,'YTick',0:5:length(nonemp));
                if z == 1
                    ylabel('Trial #');
                end
                set(gca,'fontsize',14);
                
                hold on
                
                % plot PSTH in bottom
                subplot('position',...
                    [x0+(w_stim+w_space)*(z-1) y0 w_stim h_rast]);
                
                
                psth = N_cum1{1,chann(ch)}{loc(1),loc(2),loc(3)};
                bar(T_vec,psth,'Facecolor','k','EdgeColor','k')
                xlim([-T_before T_after])
                ylim([0 Max]); box off;
                set(gca,'YTick',0:10:Max);
                if max(psth) > Max
                    ylim([0 Max+50])
                    set(gca,'YTick',0:20:(Max+50));
                end
                xlabel('Time (s)')
                if z == 1
                    ylabel('Spike Count')
                end
                set(gca,'fontsize',14);
                hold on

                %calculate FR
                tStart = 0.02;
                tEnd = 2.98;
                nTrials = length(nonemp);
                startidx = find(T_vec>=tStart,1);
                endidx = find(T_vec<=tEnd,1,'last');
                fr = sum(psth(startidx:endidx))./(nTrials*(tEnd-tStart)); %spikes per second per trial
                
                temp = cellfun(@(x) numel(x >= tStart & x < tEnd),spikes,'UniformOutput',false);
                temp = cell2mat(temp)/(tEnd-tStart);
                fr_std = std(temp);
                
                r0_t = 1.9; %total time duration before & after stimulus (exc. laser onset and offset)
                r0_fr = sum(psth([1:startidx-4 endidx+2:end]))./(nTrials*(r0_t));
                
                frStr = {['avgFR = ' num2str(fr) ' ± ' num2str(fr_std) ' Hz'];...
                         ['r_0FR = ' num2str(r0_fr) ' Hz']};
                                
                annotation('textbox',[x0+(w_stim+w_space)*(z-1) y0+h_rast-0.1 .2 .1],...
                           'string',frStr,...
                           'FitBoxToText','on',...
                           'LineStyle','none')
                
            end
            figname =sprintf('CH%g target %s%c, masker %s%c, %s',chann(ch),tdegree,char(176),mdegree,char(176),type1);
            sgtitle(figname);
            saveas(gcf,[subfolder,'/Rasters/CH',num2str(chann(ch)),'/',figname],'png');
            savefig(gcf,[subfolder,'/Rasters/CH',num2str(chann(ch)),'/',figname]);
            clf;
        end
    end
    
    % CLEAN RASTERS
    
    for i = 1:4
        for z = 1:2
                       
            loc(1) = i;   % target location
            loc(2) = z;   % target identity
            
            switch i
                case 1
                    tdegree='90';
                case 2
                    tdegree='45';
                case 3
                    tdegree='0';
                case 4
                    tdegree='-90';
            end
                        
            nonemp = find(~cellfun(@isempty,in2{1,chann(ch)}(:,i,z)));
            spikes{length(nonemp)} = [];
            
            if length(nonemp) <= 3 % for TM configurations with no data
                break
            end
            
            for x=1:length(nonemp)
                if size(in2{1,chann(ch)}{nonemp(x),loc(1),loc(2)},1) > size(in2{1,chann(ch)}{nonemp(x),loc(1),loc(2)},2)
                    spikes{x} = in2{1,chann(ch)}{nonemp(x),loc(1),loc(2)}';
                else
                    spikes{x} = in2{1,chann(ch)}{nonemp(x),loc(1),loc(2)};
                end
            end
            

            % plot stimulus on top
            subplot('position',...
                [x0+(w_stim+w_space)*(z-1) y0+(h_rast*2 +h_space*2) w_stim h_stim]);
            axis tight
            
            plot(t_stim,soundstim(:,loc(2)),'k')
            set(gca,'XTick', [],'ytick',[]);
            ylim([-1 1])
            xlim([-T_before T_after]); box off
            
            hold on
            
            % plot raster in middle
            subplot('position',...
                [x0+(w_stim+w_space)*(z-1) y0+(h_rast+h_space) w_stim h_rast]);
            
            for sp = 1:length(nonemp)
                plot([spikes{sp};spikes{sp}],[ones(size(spikes{sp}))*sp;ones(size(spikes{sp}))*sp-1],'-k','linewidth',.2)
                hold on
            end
            
            ylim([0 length(nonemp)])
            xlim([-T_before T_after])
            set(gca,'TickDir','out') % draw the tick marks on the outside
            set(gca,'XTick', []) % don't draw y-axis ticks
            %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
            %set(gca,'Color',get(gcf,'Color')) % match figure background
            %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
            box off
            set(gca,'YTick',0:5:length(nonemp));
            if z == 1
            ylabel('Trial #');
            end
            set(gca,'fontsize',14);
            
            % plot PSTH in bottom
            subplot('position',...
                [x0+(w_stim+w_space)*(z-1) y0 w_stim h_rast]);
            
            psth = N_cum2{1,chann(ch)}{loc(1),loc(2)};
            bar(T_vec,psth,'Facecolor','k','EdgeColor','k')
            xlim([-T_before T_after])
            ylim([0 Max]); box off;
            set(gca,'YTick',0:10:Max);
            if max(psth) > Max
               ylim([0 Max+50]) 
               set(gca,'YTick',0:20:(Max+50));
            end
            xlabel('Time (s)')
            if z == 1
                ylabel('Spike Count')
            end
            set(gca,'fontsize',14);
            hold on
            
            %calculate FR
            tStart = 0;
            tEnd = 2.98;
            nTrials = length(nonemp);
            startidx = find(T_vec>=tStart,1);
            endidx = find(T_vec<=tEnd,1,'last');
            fr = sum(psth(startidx:endidx))./(nTrials*(tEnd-tStart)); %spikes per second per trial

            temp = cellfun(@(x) numel(x >= tStart & x < tEnd),spikes,'UniformOutput',false);
            temp = cell2mat(temp)/(tEnd-tStart);
            fr_std = std(temp);
            
            r0_t = 1.9; %total time duration before & after stimulus (exc. laser onset and offset)
            r0_fr = sum(psth([1:startidx-4 endidx+2:end]))./(nTrials*(r0_t));

            frStr = {['avgFR = ' num2str(fr) ' ± ' num2str(fr_std) ' Hz'];...
                ['r_0FR = ' num2str(r0_fr) ' Hz']};
            
            annotation('textbox',[x0+(w_stim+w_space)*(z-1) y0+h_rast-0.1 .2 .1],...
                'string',frStr,...
                'FitBoxToText','on',...
                'LineStyle','none')
            
        end
        
        figname = sprintf('CH%g clean, %s%c, %s',chann(ch),tdegree,char(176),type1);
        sgtitle(figname);
        saveas(gcf,[subfolder,'/Rasters/CH',num2str(chann(ch)),'/',figname],'png');
        savefig(gcf,[subfolder,'/Rasters/CH',num2str(chann(ch)),'/',figname]);
        clf;
    end
end

close(gcf);

end
