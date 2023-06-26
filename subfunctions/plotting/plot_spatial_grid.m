function plot_spatial_grid(ch,nch,type,TMR,Max,max_masked,barmin,barmax)

figure;

if ~isempty(max_masked)
    %width = 4; hwratio = 1;
    width = 2.25; hwratio = 1;
else
    %width = 4; hwratio = 0.4; 
end

spatialgrid_width=0.051*4.86/2;
figuresize(width,width*hwratio,gcf,'inches')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width width*hwratio]);
colormap(parula)
xy1 = 0.25;

% Clean spots
axes('pos',[xy1 1-2*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width/hwratio])
image(Max,'CDataMapping','scaled')

caxis([barmin barmax]);

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

title(sprintf('%s %g, %gdB',type,nch,TMR),'fontweight','normal','fontsize',10);
ax = gca;
ax.FontSize = 8;

% Masked grid

if ~isempty(max_masked)
    axes('pos',[xy1 1-6.5*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width*4/hwratio])
    image(max_masked,'CDataMapping','scaled')
    
    caxis([barmin barmax]); colormap; colorbar;
    
    axis image;
    set(gca,'xtick',[1 2 3 4],'xTickLabel',{'90','45','0','-90'})
    set(gca,'ytick',[1 2 3 4],'yTickLabel',{'-90','0','45','90'})
    xlabel('Target Location (°)'); 
    ylabel('Masker Location (°)');
    
    h_bar = findobj(gcf,'Tag','Colorbar');
    set(h_bar, 'Position',[xy1+spatialgrid_width*4.25 1-6.5*spatialgrid_width/hwratio .25*spatialgrid_width spatialgrid_width*4/hwratio])
    set(h_bar,'ytick',[barmin barmax],'yticklabel',{sprintf('%g%',barmin),sprintf('%g%',barmax)},'FontSize',8);
    h_bar.Label.String = 'Performance'; h_bar.Label.Position(1) = h_bar.Label.Position(1)/2;
    h_bar.TickLabels = cellfun(@(x) [x '%'],h_bar.TickLabels,'uniformoutput',false);

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
    
    % Change text color to white for bluer spots
    
    % for i=1:4
    %     for j=1:4
    %         if mp{1,ch}((i-1)*4+j)<50
    %             text( j,i, mpstr{1,ch}{(i-1)*4+j},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
    %         end
    %     end
    % end
    
    set(gca,'fontsize',8);
end

end