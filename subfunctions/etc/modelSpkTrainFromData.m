function new_data = modelSpkTrainFromData(data)

Chans = find(~cellfun(@isempty,data.Spks_clean));

new_data = struct;

t_vec = -1:0.02:4;

for ch = Chans
    % clean trials
    for i = 1:4
        for z = 1:2
            spks = data.Spks_clean{ch}(:,i,z);
            spks{10} = [];
            new_spks = genTrains(spks);
            new_data.Spks_clean{ch}(:,i,z) = new_spks;
            new_data.n_cum_clean{ch}{i,z} = [histcounts(vertcat(new_spks{:}),t_vec),0];
        end
    end


    % masked trials
    for i = 1:4
        for j = 1:4
            for z = 1:2
                spks = data.Spks_masked{ch}(:,i,j,z);
                spks{10} = [];
                new_spks = genTrains(spks);
                new_data.Spks_masked{ch}(:,i,j,z) = new_spks;
                new_data.n_cum_masked{ch}{i,j,z} = [histcounts(vertcat(new_spks{:}),t_vec),0];
            end
        end
    end
end

new_data.Spks_clean{32} = [];
new_data.Spks_masked{32} = [];
new_data.n_cum_clean{32} = [];
new_data.n_cum_masked{32} = [];

end


function spks = genTrains(spks)

t_vec = -1:0.0025:4;

% bin non-empty trials
nonemp = find(~cellfun(@isempty,spks));

for k = 1:length(t_vec)-1
    spksPerTrial = cellfun(@(x) sum(x >= t_vec(k) & x < t_vec(k+1)),spks(nonemp));

    % calculate mean and standard deviation of spikes per bin

    mu(k) = mean(spksPerTrial);
    sigma(k) = std(spksPerTrial);
end

numNewTrials = 10 - length(nonemp);

for k = 1:length(t_vec)-1

    % generate # spikes in bin based on gaussian distribution from data

    nSpksBin = round(normrnd(mu(k),sigma(k),[1 numNewTrials]));
    nSpksBin(nSpksBin < 0) = 0;

    % assign spike times using uniform distribution within time bin
    for tt = 1:numNewTrials
        spks{length(nonemp)+tt} = [spks{length(nonemp)+tt} ; 0.0025*rand(nSpksBin(tt),1) + t_vec(k)];
    end 
end

end