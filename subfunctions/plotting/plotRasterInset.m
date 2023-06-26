function plotRasterInset(in1,N_cum1,in2,N_cum2,type1,chann,T_before,T_after,ii,jj,zz,...
    Max,subfolder)
% in1 is spike times, i.e Spks_clean{1,ch}, loc is which location i.e [speaker,song]
% from {x,z,k} = {speaker,song,trial} so c = x
% N_cum1 is the binned spike count n_cum_clean{1,ch}{speaker,song}
% (cum is summed across all trials)
% type1 is title (laser/control)
% chan is channel Max is max number of spikes in the binned
% in1,sti1,n_cum1 == MASKERS, the 2's are the clean

% plotRasterInset(in1,N_cum1,in2,N_cum2,type1,chann,T_before,T_after,...
%     ii,jj,zz,T_vec,Max,subfolder)

if ~isfolder([subfolder,'/Raster Insets'])
    mkdir([subfolder,'/Raster Insets'])
end

% call figure here for efficiency
width = 2.5; height = 4;
figure('units','inches','position',[2 2 width 4]);

if isempty(chann)
   return 
end

% define subplot heights and width
w_plot = 0.8;
h_plot = 0.38;

h_space = 0.06;

y0 = 0.15;
x0 = 0.14;

T_vec = -1:1/50:4;

for ch = 1:length(chann)
            
    disp([type1,' CH',num2str(chann(ch))]);
    
    if ~isfolder([subfolder,'/Raster Insets/CH',num2str(chann(ch))])
        mkdir([subfolder,'/Raster Insets/CH',num2str(chann(ch))])
    end
    
    % MASKED RASTERS
    
    loc = [];
    for i = ii
        for j = jj
            for z = zz
                
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
                [x0 y0+(h_plot+h_space) w_plot h_plot]);
                axis tight
                
                for sp = 1:length(nonemp)
                    plot([spikes{sp};spikes{sp}],[ones(size(spikes{sp}))*sp;ones(size(spikes{sp}))*sp-1],'-k','linewidth',.2)
                    hold on
                end
                
                ylim([0 length(nonemp)])
                xlim([-T_before T_after])
                set(gca,'TickDir','out') % draw the tick marks on the outside
                set(gca,'XTick', [],'fontsize',14) % don't draw y-axis ticks
                %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
                %set(gca,'Color',get(gcf,'Color')) % match figure background
                %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
                box off
                set(gca,'YTick',0:5:length(nonemp));
                
                % plot PSTH in bottom
                subplot('position',...
                [x0 y0 w_plot h_plot]);
                
                psth = N_cum1{1,chann(ch)}{loc(1),loc(2),loc(3)};
                plot(T_vec,psth,'k','linewidth',2)
                xlim([-T_before T_after]); set(gca,'xtick',-T_before:0.1:T_after);
                ylim([0 Max])
                set(gca,'YTick',0:10:Max,'fontsize',14);
                xlabel('Time (s)'); box off
                
                figname =sprintf('CH%g song %g T%s%c, M%s%c, %s',chann(ch),loc(3),tdegree,char(176),mdegree,char(176),type1);
                %sgtitle(figname);
                saveas(gcf,[subfolder,'/Raster Insets/CH',num2str(chann(ch)),'/',figname],'png');
                savefig(gcf,[subfolder,'/Raster Insets/CH',num2str(chann(ch)),'/',figname]);
                clf;
            end
        end
    end
    
    % CLEAN RASTERS
    
    for i = ii
        for z = zz
                       
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
            
            % plot raster on top
            subplot('position',...
                [x0 y0+(h_plot+h_space) w_plot h_plot]);
            axis tight
            
            for sp = 1:length(nonemp)
                plot([spikes{sp};spikes{sp}],[ones(size(spikes{sp}))*sp;ones(size(spikes{sp}))*sp-1],'-k','linewidth',.2)
                hold on
            end
            
            ylim([0 length(nonemp)])
            xlim([-T_before T_after])
            set(gca,'TickDir','out') % draw the tick marks on the outside
            set(gca,'XTick', [],'fontsize',14) % don't draw y-axis ticks
            %set(gca,'PlotBoxAspectRatio',[1 0.05 1]) % short and wide
            %set(gca,'Color',get(gcf,'Color')) % match figure background
            %set(gca,'YColor',get(gcf,'Color')) % hide the y axis
            box off
            set(gca,'YTick',0:5:length(nonemp));
            
            % plot PSTH in bottom
            subplot('position',...
                [x0 y0 w_plot h_plot]);
            
            psth = N_cum2{1,chann(ch)}{loc(1),loc(2)};
            plot(T_vec,psth,'k','linewidth',2)
            xlim([-T_before T_after]); set(gca,'xtick',-T_before:0.1:T_after);
            ylim([0 Max])
            set(gca,'YTick',0:10:Max,'fontsize',14);
            xlabel('Time (s)'); box off
            
            figname = sprintf('CH%g song %g, T%s%c, %s',chann(ch),loc(2),tdegree,char(176),type1);
            %sgtitle(figname);
            saveas(gcf,[subfolder,'/Raster Insets/CH',num2str(chann(ch)),'/',figname],'png');
            savefig(gcf,[subfolder,'/Raster Insets/CH',num2str(chann(ch)),'/',figname]);
            clf;

        end
        
    end
end

% close(gcf);

end
