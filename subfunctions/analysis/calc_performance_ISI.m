function [perf_file,goodchans,perf_struct] = calc_performance_ISI(Spks_clean,Spks_masked,folder,...
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
    
    data = []; groups = [];
    
    for i = 1:4 % stimuli loc
        for s = 1:2
            for k = 1:len_min{1,ch}(i)
                spTimes{1,ch}{i}{k,s} = Spks_clean{1,ch}{nonempty{1,ch}{i,s}(k),i,s}';
            end
        end
        
        if k > 1
            % SpikeTrainSet only accepts cells in row vector form, with
            % each double also in row form
            input = reshape(spTimes{1,ch}{i},1,2*k);
            STS = SpikeTrainSet(input,0,3);
            distMat = STS.ISIdistanceMatrix(0,3);
            
            [Max{1,ch}(i),~,clean_performance{1,ch}{i}] = calcpcStatic(distMat, k, 2, 0);
        else % if there's only one trial, default to chance performance
            Max{1,ch}(i) = 0;
        end
        
        
        for j = 1:4 % masker loc
            for s = 1:2 %stimuli#
                for k2 = 1:len2_min{1,ch}(i,j)
                    sp_masked{1,ch}{i,j}{k2,s} = Spks_masked{1,ch}{nonempty2{1,ch}{i,j,s}(k2),i,j,s}';
                end
            end
           
            % flip masker location for plotting purposes (j -> 5-j)
            % j = 1 (90deg) -> 4th index in 2nd dimension for
            % masked_performance = 90deg
            
            % bottom-left to top-right diagonal = co-located axis, same
            % format as in the spatial grids (since imagesc plots with
            % origin in top-left)
            
            if k2 > 1
                input = reshape(sp_masked{1,ch}{i,j},1,2*k2);
                STS = SpikeTrainSet(input,0,3);
                distMat_masked = STS.ISIdistanceMatrix(0,3);
                
                [max_masked{1,ch}(5-j,i),~,masked_performance{1,ch}{5-j,i}] = calcpcStatic(distMat_masked, k2, 2, 0);
            else
                max_masked{1,ch}(5-j,i) = 0;
            end
                        
        end
    end
    
end

if saveFile
    perf_file = [folder,'/',subject,'-',num2str(power),'mW_',num2str(TMR),...
        'dBTMR_performance_ISI',type,'.mat'];
    
    save(perf_file,'Max','max_masked','clean_performance','masked_performance','goodchans');
end

perf_struct = struct;
perf_struct.Max = Max;
perf_struct.max_masked = max_masked;

end