function [pct,ISIs] = calcISIViolations(ctrl_spks,laser_spks,ch,trialinfo)

nClean = ones(4,2);
nMasked = ones(4,4,2);
nCleanlaser = ones(4,2);
nMaskedlaser = ones(4,4,2);

TS = []; offset = 0;

for n = 1:size(trialinfo,1)
    tloc = trialinfo(n,1);
    mloc = trialinfo(n,2);
    tt = trialinfo(n,3);
    
    if trialinfo(n,4) == 0 % control
        if mloc == 0
            if nClean(tloc,tt) < 10  % clean
            temp = ctrl_spks.Spks_clean{ch}{nClean(tloc,tt),tloc,tt};
            nClean(tloc,tt) = nClean(tloc,tt) + 1;
            end
        elseif nMasked(tloc,mloc,tt) < 10 % masked
            %else
            temp = ctrl_spks.Spks_masked{ch}{nMasked(tloc,mloc,tt),tloc,mloc,tt};
            nMasked(tloc,mloc,tt) = nMasked(tloc,mloc,tt) + 1;
        end
    else % laser
        if mloc == 0 
            if nCleanlaser(tloc,tt) < 10 % clean
            temp = laser_spks.Spks_clean{ch}{nCleanlaser(tloc,tt),tloc,tt};
            nCleanlaser(tloc,tt) = nCleanlaser(tloc,tt) + 1;
            end
        elseif nMaskedlaser(tloc,mloc,tt) < 10 % masked
            %else
            temp = laser_spks.Spks_masked{ch}{nMaskedlaser(tloc,mloc,tt),tloc,mloc,tt};
            nMaskedlaser(tloc,mloc,tt) = nMaskedlaser(tloc,mloc,tt) + 1;
        end
    end
    
    TS = cat(1,TS,temp + offset);
    offset = offset + 5;
end

ISIs = diff(TS);
pct = sum(ISIs < 0.002)/(numel(TS)-1);

end