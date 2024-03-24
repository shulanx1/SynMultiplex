%% Figure 3A-C, Figure 3F-G
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
addpath(genpath('plotting'))
idx_cell = 15;  % 3; %for figure 3F&G
idx_trials = [11,8,5];  % [3,2,1]; %for figure 3F&G
dt = 5e-5;
t = -500*dt:dt:dt*(size(data(data_decode(idx_cell).idx(idx_trials(3))).EPSP_singleclust, 1)-1-500);

figure
plot(t, data(data_decode(idx_cell).idx(idx_trials(3))).EPSP_singleclust)

for i = 1:length(idx_trials)
    V = data(data_decode(idx_cell).idx(idx_trials(i))).EPSP_highcond{1};
    t = -500*dt:dt:dt*(size(V, 1)-1-500);
    figure
    subplot(211)
    for j = 1:size(V,2)
        plot(t, V(:,j)-50*(j-1), 'k')
        hold on
    end
    xlim([-0.01,0.1])
    subplot(212)
    raster_plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time(1,:))
end
%% Figure 3D
load(fullfile(pwd, 'stat', 'Figure3D.mat'))
colors = [[0, 113, 187];[0,164,79];[240,90,36]]/255;
boxplot_pairwise([datapool.pburst_1clust_lowcond', datapool.pburst_2clust_lowcond', datapool.pburst_3clust_lowcond'], colors)
ylabel(['Burst Prob.'])
boxplot_pairwise([datapool.APnum_1clust_lowcond', datapool.APnum_2clust_lowcond', datapool.APnum_3clust_lowcond'], colors)
ylabel(['AP count [100ms-1]'])
boxplot_compact({1000*datapool.isi_1clust_lowcond', 1000*datapool.isi_2clust_lowcond', 1000*datapool.isi_3clust_lowcond'}, colors)
ylabel(['ISI [ms]'])
ylim([5, 25])
%% Figure 3E-H
load(fullfile(pwd, 'stat', 'Figure3EH.mat'))
colors = [[0,164,79];[0, 113, 187]]/255;
boxplot_compact({datapool.maxarea_burst_lowcond, datapool.maxarea_single_lowcond}, colors)
boxplot_compact({datapool.sumarea_burst_lowcond, datapool.sumarea_single_lowcond}, colors)
colors = [[255, 0, 0];[0, 0, 0]]/255;
boxplot_compact({datapool.jitter_strong_lowcond, datapool.jitter_weak_lowcond}, colors)
boxplot_compact({datapool.median_onset_strong_lowcond, datapool.median_onset_weak_lowcond}, colors)
%% Figure S7
load(fullfile(pwd, 'stat', 'FigureS7.mat'))
colors = [[255, 0, 0];[255, 77, 77]; [255, 153, 153]; [0,0,0]]/255;
boxplot_with_datapoint({jitter_3strong, jitter_2strong, jitter_1strong, jitter_weak}, colors)
colors = [[153, 153, 153];[102, 102, 102]; [51, 51, 51]; [0, 0, 0]; [255, 0, 0]]/255;
boxplot_with_datapoint({jitter_1weak, jitter_2weak, jitter_3weak, jitter_4weak, jitter_4clust_strong}, colors)