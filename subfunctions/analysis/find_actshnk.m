function [min_csd,csd_onset,csd_offset,thresh_csd,gran_layers] = find_actshnk(csd,csd_t,t_resp)

% determine granular channel based on CSD sink onset latency

shnknames = {'shank1','shank2','shank3','shank4'};

allsinks = [];
allonsets = [];

t_zero = sum(csd_t < 0);

for sh = 1:4
    
    % find CSD shank from each channel + threshold for sink onset/offset
    for n = 1:6
        min_csd.(shnknames{sh})(n) = min(csd.(shnknames{sh})(n,csd_t>=t_resp(1) & csd_t<=t_resp(2)));
        thresh_csd.(shnknames{sh})(n) = -3*std(csd.(shnknames{sh})(n,csd_t < 0));
        
        csd_onset.(shnknames{sh})(n) = min([csd_t(t_zero + find(csd.(shnknames{sh})(n,csd_t >= 0) < thresh_csd.(shnknames{sh})(n),1)) Inf]);
        csd_offset.(shnknames{sh})(n) = min([csd_t(t_zero + find(csd.(shnknames{sh})(n,csd_t > csd_onset.(shnknames{sh})(n)) ...
            > thresh_csd.(shnknames{sh})(n),1)) Inf]);
        
    end
    
    % determine sink onset and offset
    
    % find channel with smallest CSD value in each shank for each tone
    [Y1.(shnknames{sh}),sinkind.(shnknames{sh})] = min(min_csd.(shnknames{sh}));
    [Y2.(shnknames{sh}),onsetind.(shnknames{sh})] = min(csd_onset.(shnknames{sh}));
    
    allsinks = cat(2,allsinks,Y1.(shnknames{sh}));
    allonsets = cat(2,allonsets,Y2.(shnknames{sh}));

    % find granular channel based on earliest onset
    sinkchan = find(csd_onset.(shnknames{sh}) == min(csd_onset.(shnknames{sh})));
    
    % if more than one channel has the minimum, determine which one to use
    % based on neighboring onsets
    neighbor = [];
    if numel(sinkchan) > 1
        for g = 1:length(sinkchan)
            if sinkchan(g) == 1
                neighbor(g) = csd_onset.(shnknames{sh})(sinkchan(g)+1)-csd_onset.(shnknames{sh})(sinkchan(g));
            elseif sinkchan(g) == 6
                neighbor(g) = csd_onset.(shnknames{sh})(sinkchan(g)-1)-csd_onset.(shnknames{sh})(sinkchan(g));
            else
                neighbor(g) = min([csd_onset.(shnknames{sh})(sinkchan(g)-1) csd_onset.(shnknames{sh})(sinkchan(g)+1)]-csd_onset.(shnknames{sh})(sinkchan(g)));
            end
        end
        [~,ind] = min(neighbor);
        gran_layers(sh) = sinkchan(ind);
    else
        gran_layers(sh) = sinkchan;
    end

    
end

end

