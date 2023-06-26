function [waveforms,IDs,groups] = sort_snippets(TS,Mua,snips,codes,labels,t_vec)

for i = 1:length(TS)
    temp = find(Mua >= (TS(i) + t_vec(1)) & Mua < (TS(i) + t_vec(end)));
    waveforms{i} = snips(temp,:);
    IDs{i} = codes(temp);
    if ~isempty(labels)
        groups{i} = labels(temp);
    end
end

end