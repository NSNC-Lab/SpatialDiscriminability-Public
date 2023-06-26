function plotSigGrids(TMR,power,folder,ctrl_perf,laser_perf,duds)

%%% INPUTS
% TMR: target-to-masker ratio
% figuresdir: subject name
% type: 'control' or 'laser'
% d1: struct of performance data

% check if faulty channels are removed
chFlag = length(ctrl_perf.d_masked) == (32 - length(duds));

nskips = 0;

for ch = 1:length(ctrl_perf.d_masked)
    
    if ~isempty(duds)
        if chFlag == 0 && ismember(ch,duds) % if dud channel isn't taken out
            continue
        elseif chFlag == 1 && ch >= duds(1)   % if dud channel is taken out
            nskips = sum(ch >= duds);
        end
    end
    
    if isempty(ctrl_perf.d{ch})
        continue
    end
        
    grid_folder = [folder,'/Comparison Grids'];
    
    if ~isfolder(grid_folder)
        mkdir(grid_folder);
    end
    
    barmin = -10;
    barmax = 10;
    
    siggrid(ch,ch+nskips,TMR,ctrl_perf,laser_perf,barmin,barmax);
        
    if exist('cluster_chs','var')
        
        h = findobj(gcf,'type','axes');  
        axes(h(2));
        title(sprintf('%s Cluster %g (CH%g), %gdB',type,cluster_ids(ch),cluster_chs(ch),TMR));

        saveas(gcf,[grid_folder '/' num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_diff.png']);
        savefig(gcf,[grid_folder '/' num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_diff.fig']);
    else
        saveas(gcf,[grid_folder '/' ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_diff.png']); 
        savefig(gcf,[grid_folder '/' ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_diff.fig']);
    end
    
    close;
end

end

function siggrid(ch,nch,TMR,ctrl_perf,laser_perf,barmin,barmax)

figure('color','w');

width = 4; hwratio = 1;

spatialgrid_width=0.051*4.86/2;
figuresize(width,width*hwratio,gcf,'inches')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width width*hwratio]);
colormap(gray)
xy1 = 0.25;

% Clean spots
axes('pos',[xy1 1-2*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width/hwratio])
image(ctrl_perf.Max{1,ch} - laser_perf.Max{1,ch},'CDataMapping','scaled')

caxis([barmin barmax]);

axis image;
ylabel('Clean','fontsize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])

for i=1:4
    mp(i)=round(ctrl_perf.Max{1,ch}(i) - laser_perf.Max{1,ch}(i));
    [~,p] = ttest2(ctrl_perf.clean_performance{1,ch}{i},laser_perf.clean_performance{1,ch}{i},...
        'vartype','unequal');
    mps{i}=makeStar(p);
end

mpstr=strtrim(cellstr(num2str(mp')));
mpsig=mps';
axis image
xg=0:1:3;
text(xg+1,[1 1 1 1], mpstr,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16)
text(xg+1,[1 1 1 1]-0.25, mpsig,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12)

for i=1:4
    if mp(i) <= 3*(barmax-barmin)/10
        text(i,1, mpstr{i},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
        text(i,1-0.25, mpsig{i},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12);
    end
end

clear mp

title(sprintf('CH%g, %gdB',nch,TMR));
ax = gca;
ax.FontSize = 14;

% Masked grid

    axes('pos',[xy1 1-6.5*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width*4/hwratio])
    image(ctrl_perf.max_masked{1,ch} - laser_perf.max_masked{1,ch},'CDataMapping','scaled')
    
    caxis([barmin barmax]); colormap; colorbar;
    
    axis image;
    set(gca,'xtick',[1 2 3 4],'xTickLabel',{'90','45','0','-90'})
    set(gca,'ytick',[1 2 3 4],'yTickLabel',{'-90','0','45','90'})
    xlabel('Target Location (°)'); 
    ylabel('Masker Location (°)');
    
    h_bar = findobj(gcf,'Tag','Colorbar');
    set(h_bar, 'Position',[xy1+spatialgrid_width*4.25 1-6.5*spatialgrid_width/hwratio .25*spatialgrid_width spatialgrid_width*4/hwratio])
    set(h_bar,'ytick',[barmin barmax],'yticklabel',{sprintf('%g%',barmin),sprintf('%g%',barmax)})
    ylabel(h_bar,'Control - Laser');
    for i=1:4 % row = masker location (1 = -90, 4 = +90)
        for j=1:4 % column = target location (1 = +90, 4 = -90)
            mp((i-1)*4+j)=round(ctrl_perf.max_masked{1,ch}(i,j) - laser_perf.max_masked{1,ch}(i,j));
            [~,p] = ttest2(ctrl_perf.masked_performance{1,ch}{i,j},laser_perf.masked_performance{1,ch}{i,j},...
                'vartype','unequal');
            mps{(i-1)*4+j}=makeStar(p);
        end
    end
    mpstr=strtrim(cellstr(num2str(mp')));
    mpsig=mps';
    axis image
    yg = 0:1:3;
    [xlbl, ylbl] = meshgrid(xg+1, yg+1);
    text( ylbl(:),xlbl(:), mpstr,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
    text( ylbl(:),xlbl(:)-0.25, mpsig,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12);
    % Change text color to white for darker spots
    
    for i=1:4
        for j=1:4
            if mp((i-1)*4+j) <= 3*(barmax-barmin)/10
                text( j,i, mpstr{(i-1)*4+j},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
                text( j,i-0.25, mpsig{(i-1)*4+j},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12);
            end
        end
    end
    
    set(gca,'fontsize',16);

end

function txt = makeStar(p)

if p < 0.001
    txt = '***';
elseif p < 0.01
    txt = '**';
elseif p < 0.05
    txt = '*';
else
    txt = '';
end

end