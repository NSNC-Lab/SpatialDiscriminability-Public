function dist_mat = calcSpkCtDist(spks)

% spks: cell vector with each cell containing a spike train 

nT = numel(spks);

dist_mat = zeros(nT,nT);

for n1 = 1:nT
    for n2 = n1:nT
        dist_mat(n1,n2) = abs(numel(spks{n1})-numel(spks{n2}));
        dist_mat(n2,n1) = dist_mat(n1,n2);
    end
end

end