function [inds,d_FR] = findColdspots(perfs,d_vals,spk_struct,perf_thresh,d_thresh)

% subject_struct: exps(ns).ctrl_perf_SPIKE.Max or
% exps(ns).ctrl_perf_SPIKE.max_masked

% inds_final must be a column  vector

% first, find spots that fulfill both performance and d' thresholds

for n = 1:numel(perfs)
    % inds = find(perfs < perf_thresh & d_vals < d_thresh);
    inds = find(perfs < perf_thresh);
end

end