% run this loop to extract clusters to Matlab, then run analysis per subject

% subjects = {'616283','616299','616302','616303','616304','616305','616306','616307','616312'};
subjects = {'616283'};

whichpower = 10;
whichTMR = 0;
targetdB = 75;

addpath('stimuli');
addpath(genpath('subfunctions'));

nExps = length(subjects);

for ns = 1:length(subjects)

    % analysis_storage - where to store the data sorted by trial and analysis results
        
    % data_storage - where to store the spikes from kilosort

    analysis_storage = 'Arch-data/analysis'; %uigetdir(cd,'Select where to save processed spikes and results.');   % 
    data_storage = 'Arch-data/processed-Kilosort'; %uigetdir(cd,'Select folder with stored data tanks.');     
    
    addpath(analysis_storage);

    %% Loop through processed clusters to get spikes / discriminability

    disp(['Running experiment ' num2str(ns) ' out of ' num2str(nExps) '.']);

    TMR = whichTMR;
    power = whichpower;
    subject = subjects{ns};

    ksort_folder = fullfile(data_storage,subject,'Kilosort-merged');

    spgrid_Phy_to_results;

    clearvars ctrl_spks ctrl_perf laser_spks laser_perf
    close all;
    clc;

end

% next: finalPopAnalysis
