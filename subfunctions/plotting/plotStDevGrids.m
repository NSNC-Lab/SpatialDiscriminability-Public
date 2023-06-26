function plotStDevGrids(TMR,power,subfolder,type,d1,duds,cluster_chs,cluster_ids)

%%% INPUTS
% TMR: target-to-masker ratio
% figuresdir: subject name
% type: 'control' or 'laser'
% d1: struct of performance data

% check if faulty channels are removed
chFlag = length(d1.d_masked) == (32 - length(duds));

nskips = 0;

for ch = 1:length(d1.d_masked)
    
    if ~isempty(duds)
        if chFlag == 0 && ismember(ch,duds) % if dud channel isn't taken out
            continue
        elseif chFlag == 1 && ch >= duds(1)   % if dud channel is taken out
            nskips = sum(ch >= duds);
        end
    end
    
    if isempty(d1.d{ch})
        continue
    end
        
    grid_folder = [subfolder,'/StDev Grids'];
    
    if ~isfolder(grid_folder)
        mkdir(grid_folder);
    end
    
    barmin = 0;
    barmax = 3;
    
    stdevgrid(ch+nskips,type,TMR,d1.d{1,ch},d1.d_masked{1,ch},...
        d1.p{1,ch},d1.p_masked{1,ch},barmin,barmax)
        
    if exist('cluster_chs','var')
        
        h = findobj(gcf,'type','axes');  
        axes(h(2));
        title(sprintf('%s Cluster %g (CH%g), %gdB',type,cluster_ids(ch),cluster_chs(ch),TMR));

        saveas(gcf,[grid_folder '/' type num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE.png']);
%         savefig(gcf,[grid_folder '/' type num2str(cluster_ids(ch)) '(CH' num2str(cluster_chs(ch))  '), ' num2str(TMR) 'dBTMR, ' ...
%             num2str(power) 'mW_SPIKE.fig']);
    else
        saveas(gcf,[grid_folder '/' type ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
            num2str(power) 'mW_SPIKE.png']); 
%         savefig(gcf,[grid_folder '/' type ', CH' num2str(ch+nskips) ', ' num2str(TMR) 'dBTMR, ' ...
%             num2str(power) 'mW_SPIKE.fig']);
    end
    
    close;
end

end

function stdevgrid(nch,type,TMR,d,d_masked,p,p_masked,barmin,barmax)

figure('color','w');
spatialgrid_width=0.051*4.86/2;
figuresize(width,width*hwratio,gcf,'inches')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width width*hwratio]);

if ~isempty(d_masked)
    width = 4; hwratio = 1;
else
    width = 4; hwratio = 0.4; 
end

spatialgrid_width=0.051*4.86/2;
figuresize(width,width*hwratio,gcf,'inches')
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 width width*hwratio]);
colormap(gray)
xy1 = 0.25;

% Clean spots
axes('pos',[xy1 1-2*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width/hwratio])
image(abs(d),'CDataMapping','scaled')

caxis([barmin barmax]);

axis image;
ylabel('Clean','fontsize',16)
set(gca,'xtick',[])
set(gca,'ytick',[])

for i=1:4
    mp(i)=round(10*d(i))/10;
    mps{i}=makeStar(p(i));
end

mpstr=strtrim(cellstr(num2str(mp')));
mpsig=mps';
axis image
xg=0:1:3;
text(xg+1,[1 1 1 1], mpstr,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16)
text(xg+1,[1 1 1 1]-0.25, mpsig,'color','k','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12)

for i=1:4
    if abs(mp(i)) <= 3*barmax/10
        text(i,1, mpstr{i},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
        text(i,1-0.25, mpsig{i},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12);
    end
end

clear mp

title(sprintf('%s %g, %gdB',type,nch,TMR));
ax = gca;
ax.FontSize = 14;

% Masked grid

if ~isempty(d_masked)
    axes('pos',[xy1 1-6.5*spatialgrid_width/hwratio spatialgrid_width*4 spatialgrid_width*4/hwratio])
    image(abs(d_masked),'CDataMapping','scaled')
    
    caxis([barmin barmax]); colormap; colorbar;
    
    axis image;
    set(gca,'xtick',[1 2 3 4],'xTickLabel',{'90','45','0','-90'})
    set(gca,'ytick',[1 2 3 4],'yTickLabel',{'-90','0','45','90'})
    xlabel('Target Location (°)'); 
    ylabel('Masker Location (°)');
    
    h_bar = findobj(gcf,'Tag','Colorbar');
    set(h_bar, 'Position',[xy1+spatialgrid_width*4.25 1-6.5*spatialgrid_width/hwratio .25*spatialgrid_width spatialgrid_width*4/hwratio])
    set(h_bar,'ytick',[barmin barmax],'yticklabel',{sprintf('%g%',barmin),sprintf('%g%',barmax)})
    for i=1:4
        for j=1:4
            mp((i-1)*4+j)=round(10*d_masked(i,j))/10;
            mps{(i-1)*4+j}=makeStar(p_masked(i,j));
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
            if abs(mp((i-1)*4+j)) <= 3*barmax/10
                text( j,i, mpstr{(i-1)*4+j},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',16);
                text( j,i-0.25, mpsig{(i-1)*4+j},'color','white','HorizontalAlignment','center','VerticalAlignment','middle','fontsize',12);
            end
        end
    end
    
    set(gca,'fontsize',16);
end

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