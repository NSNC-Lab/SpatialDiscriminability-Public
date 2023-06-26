function plotGrids(ctrl_perf,laser_perf,TMR,power,subfolder)

% Plots control and laser grids together

for ch = 1:length(ctrl_perf.Max)
        
    if isempty(ctrl_perf.Max{ch})
        continue
    end
        
    grid_folder = [subfolder,'/Control vs Laser Grids'];
    
    if ~isfolder(grid_folder)
        mkdir(grid_folder);
    end
    
    barmin = 50; barmax = 90;
    
    width = 2.25; hwratio = 1;
    figuresize(2*width,width*hwratio,gcf,'inches')

    % control
    plotSingleGrid('Control',ctrl_perf.Max{1,ch},ctrl_perf.max_masked{1,ch},barmin,barmax,hwratio,0.05)
        
    % laser
    plotSingleGrid('Laser',laser_perf.Max{1,ch},laser_perf.max_masked{1,ch},barmin,barmax,hwratio,0.5)

    sgtitle(['CH' num2str(ch)]);

    saveas(gcf,[grid_folder '/CH' num2str(ch) ' ' num2str(TMR) 'dBTMR, ' ...
        num2str(power) 'mW.png']);
    savefig(gcf,[grid_folder '/CH' num2str(ch) ' ' num2str(TMR) 'dBTMR, ' ...
        num2str(power) 'mW']);
    close;
end

end

function plotSingleGrid(type,Max,max_masked,barmin,barmax,hwratio,xy1)

spatialgrid_width=0.051*4.86/2;
set(gcf, 'Unit', 'inches');
%set(gcf, 'PaperPosition', [0 0 width width*hwratio]);

% Clean spots
axes('pos',[xy1 1-2*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width/hwratio])
image(Max,'CDataMapping','scaled')

clim([barmin barmax]); colormap(parula)

axis image;
ylabel('Clean','fontsize',8);
set(gca,'xtick',[])
set(gca,'ytick',[])

for i=1:4
    mp(i)=round(Max(i));
end

mpstr=strtrim(cellstr(num2str(mp')));
axis image
xg=0:1:3;
text(xg+1,[1 1 1 1], mpstr,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',10)
clear mp

title(sprintf('%s',type),'fontweight','normal','fontsize',10);
ax = gca;
ax.FontSize = 8;

% Masked grid

if ~isempty(max_masked)
    axes('pos',[xy1 1-6.5*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width*4/hwratio])
    image(max_masked,'CDataMapping','scaled')
    
    clim([barmin barmax]); colormap(parula); 
        if strcmp(type,'Laser')
        colorbar; end

    axis image;
    set(gca,'xtick',[1 2 3 4],'xTickLabel',{'90','45','0','-90'})
    set(gca,'ytick',[1 2 3 4],'yTickLabel',{'-90','0','45','90'})
    xlabel('Target Location (°)'); 
    ylabel('Masker Location (°)');

    for i=1:4
        for j=1:4
            mp((i-1)*4+j)=round(max_masked(i,j));
        end
    end

    mpstr=strtrim(cellstr(num2str(mp')));
    axis image
    yg = 0:1:3;
    [xlbl, ylbl] = meshgrid(xg+1, yg+1);
    text( ylbl(:),xlbl(:), mpstr,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',10);

    set(gca,'fontsize',8);

    if strcmp(type,'Laser')
        colorbar;
        h_bar = findobj(gcf,'Tag','Colorbar');
        set(h_bar, 'Position',[xy1+spatialgrid_width*4.25 1-6.5*spatialgrid_width/hwratio .25*spatialgrid_width spatialgrid_width*4/hwratio])
        set(h_bar,'ytick',[barmin barmax],'yticklabel',{sprintf('%g%',barmin),sprintf('%g%',barmax)},'FontSize',8);
        h_bar.Label.String = 'Performance'; h_bar.Label.Position(1) = h_bar.Label.Position(1)/2;
        h_bar.TickLabels = cellfun(@(x) [x '%'],h_bar.TickLabels,'uniformoutput',false);
    end

end

end