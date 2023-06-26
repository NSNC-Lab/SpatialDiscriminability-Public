function [perf_file,goodchans,perf_struct] = calc_performance_spkct(Spks_clean,Spks_masked,folder,...
    subject,type,duds,TMR,power,saveFile)

if ~exist('saveFile','var')
   saveFile = 1; 
end

chlim = length(Spks_clean);

% determine minimum # of trials between target songs for each grid
% spot/channel
for ch = 1:chlim
    
    if chlim == 32 && ismember(ch,duds) || isempty(Spks_clean{1,ch})
        for i = 1:4 % stimuli loc
            for j = 1:4 % masker loc
                len_min{1,ch}(i) = 0;
                len2_min{1,ch}(i,j) = 0;
            end
        end
        continue
    end
    
    for i = 1:4 % stimuli loc
        for j = 1:4 % masker loc
            for s = 1:2 % stimuli#
                nonempty{1,ch}{i,s}=find(~cellfun(@isempty,Spks_clean{1,ch}(:,i,s)));
                len{1,ch}{i,s}=length(nonempty{1,ch}{i,s});
                nonempty2{1,ch}{i,j,s}=find(~cellfun(@isempty,Spks_masked{1,ch}(:,i,j,s)));
                len2{1,ch}{i,j,s}=length(nonempty2{1,ch}{i,j,s});
            end
            len_min{1,ch}(i)=min(len{1,ch}{i,:});
            len2_min{1,ch}(i,j)=min(len2{1,ch}{i,j,:});
        end
    end
    
end

goodchans = [];

for ch = 1:chlim
    
    if mean([len_min{1,ch}(:);len2_min{1,ch}(:)]) < 5 || isnan(mean([len_min{1,ch}(:);len2_min{1,ch}(:)]))
        continue
    end
    
    disp(['Calculating performance for channel ',num2str(ch)]);
    
    for i = 1:4 % stimuli loc
        for j = 1:4 % masker loc
            for s = 1:2 %stimuli#
                for k = 1:len_min{1,ch}(i)
                    temp = Spks_clean{1,ch}{nonempty{1,ch}{i,s}(k),i,s};
                    spTimes{1,ch}{i}{k,s} = temp(temp >= 0 & temp < 3);
                end
                for k2 = 1:len2_min{1,ch}(i,j)
                    temp2 = Spks_masked{1,ch}{nonempty2{1,ch}{i,j,s}(k2),i,j,s};
                    sp_masked{1,ch}{i,j}{k2,s} = temp2(temp2 >= 0 & temp2 < 3);
                end
            end
            
            
            if k > 1
                distMat = calcSpkCtDist(spTimes{1,ch}{i});
                Max{1,ch}(i) = calcpcStatic(distMat, k, 2, 0);
            else % if there's only one trial, default to chance performance
                Max{1,ch}(i) = 0;
            end
            
            % flip masker location for plotting purposes (j -> 5-j)
            % j = 1 (90deg) -> 4th index in 2nd dimension for
            % masked_performance = 90deg
            
            % bottom-left to top-right diagonal = co-located axis, same
            % format as in the spatial grids (since imagesc plots with
            % origin in top-left)
            
            if k2 > 1
                distMat_masked = calcSpkCtDist(sp_masked{1,ch}{i,j});
                max_masked{1,ch}(5-j,i) = calcpcStatic(distMat_masked,k2, 2, 0);
            else
                max_masked{1,ch}(5-j,i) = 0;
            end
            
        end
    end
        
end

if saveFile
    perf_file = [folder,'/',subject,'-',num2str(power),'mW_',num2str(TMR),...
        'dBTMR_performance_spkct',type,'.mat'];
    
    save(perf_file,'Max','max_masked');
end

perf_struct = struct;
perf_struct.Max = Max;
perf_struct.max_masked = max_masked;

end