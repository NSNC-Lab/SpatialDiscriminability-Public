function [raster,n_cum] = cal_raster(TS,Mua,t_vec)

dummy = [];
for i = 1:length(TS)
    %raster{i} = Mua(Mua >= TS(i) + t_vec(1) & Mua < TS(i) + t_vec(end)) - TS(i);
    raster{i} = sort(Mua(Mua >= TS(i) + t_vec(1) & Mua < TS(i) + t_vec(end)) - TS(i),'ascend');
    dummy = cat(1,dummy,raster{i});
end

n_cum = histcounts(dummy,t_vec);
n_cum(end+1) = 0;

end