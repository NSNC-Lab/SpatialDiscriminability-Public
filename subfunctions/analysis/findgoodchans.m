function chans = findgoodchans(perf_struct,perf_thresh,d_thresh)

% channels must have at least one spot that is above performance and
% d' thresholds

chans = find( (cell2mat(cellfun(@(x) any(max(x) > perf_thresh),perf_struct.Max,'UniformOutput',false)) &...
    cell2mat(cellfun(@(x) any(max(x) > d_thresh),perf_struct.d,'UniformOutput',false)) ) | ...
    (cell2mat(cellfun(@(x) any(max(x) > perf_thresh),perf_struct.max_masked,'UniformOutput',false)) & ...
    (cell2mat(cellfun(@(x) any(max(x) > d_thresh),perf_struct.d_masked,'UniformOutput',false)))) );

end