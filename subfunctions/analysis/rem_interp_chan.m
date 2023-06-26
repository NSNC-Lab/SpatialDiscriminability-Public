function [avg_interp] = rem_interp_chan(avg,bad_chan_ind)

avg_interp=avg;
if bad_chan_ind > 1 & bad_chan_ind < 8
    
    good_chans = [avg{bad_chan_ind-1};avg{bad_chan_ind+1}];
    
    [X,Y] = meshgrid(1:length(avg{bad_chan_ind}),1:2:3);
    [XI,YI] = meshgrid(1:length(avg{bad_chan_ind}),1:1:3);
    chan_interp=interp2(X,Y,good_chans,XI,YI,'linear');
    avg_interp{bad_chan_ind}=chan_interp(2,:);
    
end


if bad_chan_ind == 1
    avg_interp{bad_chan_ind}=avg{bad_chan_ind+1};
end

if bad_chan_ind > 7
    avg_interp{bad_chan_ind}=avg{bad_chan_ind-1};
end

end

