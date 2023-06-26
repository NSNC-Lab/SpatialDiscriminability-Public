function [allTS_new,Mua_new,snips_new,codes_new,lfp_new] = stitchSpkData(allTS,Mua,snips,codes,explens)

% combines spiking data across blocks

nBlocks = length(allTS);

[allTS_new,Mua_new,snips_new,codes_new] = ...
    concatenate_two_blocks(allTS,Mua,snips,codes,explens(1:end-1));

if nBlocks > 2
    for nT = 3:nBlocks
        
        allTS_temp = [allTS_new,allTS{nT}];
        Mua_temp = [Mua_new,allTS{nT}];
        snips_temp = [snips_new,allTS{nT}];
        codes_temp = [codes_new,allTS{nT}];
        
        [allTS_new,Mua_new,snips_new,codes_new] = ...
            concatenate_two_blocks(allTS_temp,Mua_temp,snips_temp,codes_temp,sum(explens(1:end-nT+1)));
        
    end
end

end

function [allTS_new,Mua_new,snips_new,codes_new] = ...
    concatenate_two_blocks(allTS,Mua,snips,codes,t_adjust)

nChans = length(Mua{1});

TSfields = fieldnames(allTS{1});

for ts = 1:length(TSfields)
    for n = 1:numel(allTS{1}.(TSfields{ts}))
        [i1,i2,i3] = ind2sub(size(allTS{1}.(TSfields{ts})),n);
        allTS_new.(TSfields{ts}){i1,i2,i3} = cat(1,allTS{1}.(TSfields{ts}){n},allTS{2}.(TSfields{ts}){n} + t_adjust);
    end
end

for ch = 1:nChans
    
    Mua_new{ch} = cat(1,Mua{1}{ch},Mua{2}{ch} + t_adjust);
    snips_new{ch} = cat(1,snips{1}{ch},snips{2}{ch});
    codes_new{ch} = cat(1,codes{1}{ch},codes{2}{ch});

end
    
end
