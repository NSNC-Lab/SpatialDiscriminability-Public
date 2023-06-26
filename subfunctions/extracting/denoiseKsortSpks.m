function [Mua,snips,codes,labels] = denoiseKsortSpks(Mua_in,snips_in,codes_in,labels_in)

Mua_chans = find(~cellfun(@isempty,Mua_in));

% this function denoises spikes by cluster in each channel

clen = numel(Mua_in);

Mua = cell(1,clen);
snips = cell(1,clen);
codes = cell(1,clen);
labels = cell(1,clen);

for ch = Mua_chans
    
    % remove spikes from clusters if they are artifact (absolute value of
    % above 1500microV or positive peak above 750microV)
    bad_inds = find((any(abs(snips_in{ch}) > 1500E-06,2) | any(snips_in{ch} > 750E-06,2)) | isnan(Mua_in{ch}));
    
    Mua_in{ch}(bad_inds) = [];
    snips_in{ch}(bad_inds,:) = [];
    codes_in{ch}(bad_inds) = [];
    labels_in{ch}(bad_inds) = [];
    
    Mua{ch} = Mua_in{ch};
    snips{ch} = snips_in{ch};
    codes{ch} = codes_in{ch};
    labels{ch} = labels_in{ch};
    
end
    
end