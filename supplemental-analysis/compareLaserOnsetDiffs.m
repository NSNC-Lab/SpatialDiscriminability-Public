% For Supplementary Figure S5, after getting onset PSTHs for NS cells in
% aggregateSubData and aggregateControlData

arch = load('Arch_subject_NS_onset.mat');
nonarch = load('Control_subject_NS_onset.mat');

onset = find(arch.T_vec == -0.05);

arch_diff = arch.psth_ctrl(:,onset) - arch.psth_laser(:,onset);
nonarch_diff = nonarch.psth_ctrl(:,onset) - nonarch.psth_laser(:,onset);

[~,p] = ttest2(arch_diff,nonarch_diff)
d = computeCohen_d(nonarch_diff, arch_diff,'independent');