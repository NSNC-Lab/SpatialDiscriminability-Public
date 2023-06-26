function [perf_file,goodchans,perf_struct] = calc_performance_SPIKE(Spks_clean,Spks_masked,folder,...
    subject,type,duds,TMR,power,saveFile)

% outputs:
% perf_file - matfile with performance results
% goodchans - channels with hotspots
% perf_struct - struct of performance results (saved in perf_file, but also
% called here as an output so that you don't have to load it after this calling this function)

if ~exist('saveFile','var')
    saveFile = 1;
end

perf_file = [folder,'/',subject,'-',num2str(power),'mW_',num2str(TMR),...
    'dBTMR_performance_SPIKE',type,'.mat'];

% if isfile(perf_file)
% 
%     load(perf_file);
% 
%     perf_struct = struct;
%     perf_struct.Max = Max;
%     perf_struct.max_masked = max_masked;
%     perf_struct.d = d;
%     perf_struct.d_masked = d_masked;
%     perf_struct.p = p;
%     perf_struct.p_masked = p_masked;
%     perf_struct.clean_performance = clean_performance;
%     perf_struct.masked_performance = masked_performance;
% 
%     return
% 
% end

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

locs = [90 45 0 -90];

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
            distMat = STS.SPIKEdistanceMatrix(0,3);
            % [Max{1,ch}(i),~,~,clean_performance] = calcpc(distMat, k, 2, 1,[], 'new');
            [Max{1,ch}(i),~,clean_performance{1,ch}{i}] = calcpcStatic(distMat, k, 2, 0);  % despite being called max it's actually the mean

            null_clean{1,ch}{i} = [];
            for s = 1:2
                input = repmat(spTimes{1,ch}{i}(:,s)',1,2);
                STS = SpikeTrainSet(input,0,3);
                distMat = STS.SPIKEdistanceMatrix(0,3);
                % [~,~,~,temp] = calcpcSelf(distMat, k, 2, 1,[], 'new');
                [~,~,temp] = calcpcStatic(distMat, k, 2, 1);

                null_clean{1,ch}{i} = cat(1,null_clean{1,ch}{i},squeeze(temp));
            end
            d{1,ch}(i) = computeCohen_d(clean_performance{1,ch}{i}, null_clean{1,ch}{i}, 'independent');
            data = cat(1,data,clean_performance{1,ch}{i});
            groups = cat(1,groups,repmat({['C' num2str(locs(i))]},numel(clean_performance{1,ch}{i}),1));

        else % if there's only one trial, default to chance performance
            Max{1,ch}(i) = 0;
            d{1,ch}(i) = 0;
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
                distMat_masked = STS.SPIKEdistanceMatrix(0,3);

                % [max_masked{1,ch}(5-j,i),~,~,masked_performance] = calcpc(distMat_masked,k2, 2, 1,[], 'new');
                [max_masked{1,ch}(5-j,i),~,masked_performance{1,ch}{5-j,i}] = calcpcStatic(distMat_masked, k2, 2, 0);

                null_masked{1,ch}{5-j,i} = [];
                for s = 1:2
                    input = repmat(sp_masked{1,ch}{i,j}(:,s)',1,2);
                    STS = SpikeTrainSet(input,0,3);
                    distMat_masked = STS.SPIKEdistanceMatrix(0,3);
                    % [~,~,~,temp] = calcpcSelf(distMat_masked, k2, 2, 1,[], 'new');
                    [~,~,temp] = calcpcStatic(distMat_masked, k2, 2, 1);

                    null_masked{1,ch}{5-j,i} = cat(1,null_masked{1,ch}{5-j,i},squeeze(temp));
                end
                d_masked{1,ch}(5-j,i) = computeCohen_d(masked_performance{1,ch}{5-j,i}, null_masked{1,ch}{5-j,i}, 'independent');
                data = cat(1,data,masked_performance{1,ch}{5-j,i});
                groups = cat(1,groups,repmat({['T' num2str(locs(i)) 'M' num2str(locs(j))]},numel(masked_performance{1,ch}{5-j,i}),1));

            else
                max_masked{1,ch}(5-j,i) = 0;
                d_masked{1,ch}(5-j,i) = 0;
            end


            if i == j
                ref_pt(i) = max_masked{1,ch}(5-j,i);
            end

        end
    end

    [p_anova,~,stats] = anova1(data,groups,'off');
    gnames = stats.gnames;
    if p_anova < 0.05
        temp = multcompare(stats);
        % use point closest to mean co-located performance as reference pt
        % for post-hoc comparisons
        [~,ind] = min(ref_pt - mean(ref_pt));
        ref_ind = 2+6*(ind-1);
        % keep comparisons with reference group
        post_hoc = temp(any(temp(:,[1 2]) == ref_ind,2),:);

        all_grps = setxor(1:20,ref_ind);

        p{1,ch} = ones(1,4);
        p_masked{1,ch} = ones(4,4);
        for grp = 1:size(post_hoc,1)
            [i,j,maskedFlag] = detPosition(all_grps(grp));
            if ~maskedFlag
                p{1,ch}(i) = post_hoc(grp,end);
            else
                p_masked{1,ch}(5-j,i) = post_hoc(grp,end);
            end
        end
        % organize p-values from post-hoc comparisons by clean and masked
    else
        p{1,ch} = ones(1,4);
        p_masked{1,ch} = ones(4,4);
    end

    % if channel has hotspot in both clean and masked trials, add to good
    % channels
    if any(Max{1,ch} >= 70) || any(max_masked{1,ch}(:) >= 70)
        goodchans = cat(2,goodchans,ch);
    end

end

if saveFile
    save(perf_file,'Max','max_masked','p','p_masked','d','d_masked',...
        'clean_performance','masked_performance','goodchans','null_clean','null_masked');
end

perf_struct = struct;
perf_struct.Max = Max;
perf_struct.max_masked = max_masked;
perf_struct.d = d;
perf_struct.d_masked = d_masked;
perf_struct.p = p;
perf_struct.p_masked = p_masked;
perf_struct.clean_performance = clean_performance;
perf_struct.masked_performance = masked_performance;

end

function [tloc,mloc,maskedFlag] = detPosition(grp_ind)

if rem(grp_ind-1,5) == 0 % clean
    tloc = 1 + (grp_ind-1)/5;
    mloc = 0;
    maskedFlag = 0;
else % masked
    [mloc,tloc] = ind2sub([4 4],grp_ind - ceil(grp_ind/5));
    maskedFlag = 1;
end

end
