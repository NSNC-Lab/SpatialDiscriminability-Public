function plotPopulationEnv_V2(exps,exc_cells,TMR,power,folder)

ct = 1;
for ns = 1:length(exps)
    
    %chans_ctrl = findgoodchans(exps(ns).ctrl_perf_SPIKE,perf_thresh,d_thresh);
    %chans_laser = findgoodchans(exps(ns).laser_perf_SPIKE,perf_thresh,d_thresh);
    %Chans = union(chans_ctrl,chans_laser);
    Chans = exc_cells{ns};
    
    if isempty(Chans), continue; end
    
    for ch = Chans
                
        for i = 1:4
            for j = 1:4
                mp_ctrl{1,ct}((j-1)*4+(5-i)) = exps(ns).ctrl_perf_SPIKE.max_masked{ch}(i,j);
                mp_laser{1,ct}((j-1)*4+(5-i)) = exps(ns).laser_perf_SPIKE.max_masked{ch}(i,j);
            end
        end
        ct = ct + 1;
    end
    
end
    
locs = fliplr([-90 0 45 90]);

figure('unit','inches','position',[6 3 3 2.5]); hold on;

loc_max_ctrl{4,4} = [];
loc_max_laser{4,4} = [];

% per unit
for t = 1:ct-1
    
    for i=1:4
        for j=1:4
            
            loc_max_ctrl{i,j} = [loc_max_ctrl{i,j} mp_ctrl{1,t}((j-1)*4+(5-i))];
            loc_max_laser{i,j} = [loc_max_laser{i,j} mp_laser{1,t}((j-1)*4+(5-i))];
            
            if mp_ctrl{1,t}(i+(j-1)*4) >= 70
                scatter3(locs(j),locs(i),mp_ctrl{1,t}(i+(j-1)*4),16,[0.8 0.8 0.8],'filled');
            else
                scatter3(locs(j),locs(i),mp_ctrl{1,t}(i+(j-1)*4),16,[0.8 0.8 0.8]);
            end
            
        end
    end
end

% finalxls = [header, num2cell(tempxls)];
% xlswrite(['Population results/' type ' hotspot units performance envelope, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR'],tempxls);

% determine max perfs at each site for upper envelope
loc_max2_ctrl = zeros(4,4);
loc_max2_laser = zeros(4,4);
for i=1:4
    for j=1:4
        [loc_max2_ctrl(i,j)] = max(loc_max_ctrl{i,j});
        [loc_max2_laser(i,j)] = max(loc_max_laser{i,j});

%         if loc_max2_ctrl(i,j) >= 70
%         scatter3(locs(j),locs(i),loc_max2_ctrl(i,j),16,'k','filled');
%         end
    end
end

cleanaxis = repmat(locs,4,1); %(4:-1:1)'*ones(1,4);
maskedaxis = repmat(fliplr(locs)',1,4);

C0 = zeros(4,4,3); C0(:,:,1) = 1; C0(:,:,2) = 1;

% surf(cleanaxis,maskedaxis,loc_max2_ctrl,C0,'FaceAlpha',0.3,'FaceColor','interp'); hold on;
% surf(cleanaxis,maskedaxis,loc_max2_laser,C1,'FaceAlpha',0.3,'FaceColor','interp');

loc_max2_ctrl(loc_max2_ctrl < 70) = 70;

% control upper envelope surface
surf(cleanaxis,maskedaxis,loc_max2_ctrl,C0,'FaceAlpha',0.3,'FaceColor','interp'); hold on;

% threshold surface
Ct = zeros(4,4,3); Ct(:,:,3) = 1;
surf(cleanaxis,maskedaxis,repmat(70,4,4),Ct,'FaceAlpha',0.18,'FaceColor','interp'); 

grid on
set(gca,'view',[-42 20],'ztick',[60:10:90],'gridlinestyle','--','gridalpha',0.8,...
    'xtick',[-90 0 45 90],'ytick',[-90 0 45 90],'xdir','reverse','ydir','reverse');
xtickformat('%gº'); ytickformat('%gº'); ztickformat('percentage');

zlim([60 90]); pbaspect([1 1 1.3]);

xlabel('Target location'); ylabel('Masker location'); zlabel('Performance');

% title('Upper envelope of performance','Fontsize',10);
set(gca,'FontSize',8);

savefig(gcf,[folder filesep 'Control performance envelope, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR.fig']);
saveas(gcf,[folder filesep 'Control performance envelope, ' num2str(power) 'mW ' num2str(TMR) 'dB TMR.png']);

end