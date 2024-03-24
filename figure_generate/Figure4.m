%% Figure 4A
load(fullfile(pwd, 'stat', 'Figure4A.mat'))
addpath(genpath(fullfile(pwd, 'plotting')))
colors = [[0, 0, 0]; [1, 0, 0]];
boxplot_pairwise(Rin, colors);
%% Figure 4C
load(fullfile(pwd, 'stat', 'Figure4C.mat'))
ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
curves{1} = fit([2:8]', mean(spike_prop_low)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, 1, 1, 5], 'Lower', [2, -inf, 0, 0.1],'Exclude', [7]);
curves{2} = fit([2:8]', mean(spike_prop_high)', ft,'Start', [4, 1, 0,1], 'Upper', [8, 1, 1, 5], 'Lower', [2, -inf, 0, 2],'Exclude', [4, 7]);
colors = [[0, 0, 0]; [1, 0, 0]];
errorbar_with_fitcurve([2:8], {spike_prop_low', spike_prop_high'}, curves, colors)
xlim([2,7.2])
%% Figure 4D
load(fullfile(pwd, 'stat', 'Figure4D.mat'))
colors = [[1, 0, 0]; [0, 0, 0]];
boxplot_compact({datapool.jitter_1clust_strong, datapool.jitter_1clust_weak}, colors)
%% Figure 4E
load(fullfile(pwd, 'stat', 'Figure4E.mat'))
colors = [[0, 0, 0]; [1, 0, 0]];
boxplot_pairwise(onset_pool(:,1:2), colors)
%% Figure 4G-H
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
idx_cell = 6;  
idx_trials = [4,2,1];  
dt = 5e-5;
t = -500*dt:dt:dt*(size(data(data_decode(idx_cell).idx(idx_trials(3))).EPSP_singleclust, 1)-1-500);

figure
plot(t, data(data_decode(idx_cell).idx(idx_trials(3))).EPSP_singleclust)

% low cond
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

% high cond
for i = 1:length(idx_trials)
    V = data(data_decode(idx_cell).idx(idx_trials(i))).EPSP_highcond{2};
    t = -500*dt:dt:dt*(size(V, 1)-1-500);
    figure
    subplot(211)
    for j = 1:size(V,2)
        plot(t, V(:,j)-50*(j-1), 'k')
        hold on
    end
    xlim([-0.01,0.1])
    subplot(212)
    raster_plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time(2,:))
end
%% Figure 4I
load(fullfile(pwd, 'stat', 'Figure4I.mat'))
colors = [[255, 0, 0];[0, 0, 0]]/255;
boxplot_compact({datapool.jitter_strong_highcond, datapool.jitter_weak_highcond}, colors)
boxplot_compact({datapool.median_onset_strong_highcond, datapool.median_onset_weak_highcond}, colors), ylim([0,30])
%% Figure 4J
load(fullfile(pwd, 'stat', 'Figure4J.mat'))
colors = [[0,164,79];[0, 113, 187]]/255;
boxplot_compact({datapool.maxarea_burst_highcond_highI, datapool.maxarea_single_highcond_highI}, colors)
boxplot_compact({datapool.sumarea_burst_highcond_highI, datapool.sumarea_single_highcond_highI}, colors)
%% Figure 4K
load(fullfile(pwd, 'stat', 'Figure3D.mat'))
colors = [[0, 113, 187];[0,164,79];[240,90,36]]/255;
boxplot_compact([datapool.pburst_1clust_lowcond', datapool.pburst_2clust_lowcond', datapool.pburst_3clust_lowcond'], colors)
ylabel(['Burst Prob.'])
boxplot_compact([datapool.APnum_1clust_lowcond', datapool.APnum_2clust_lowcond', datapool.APnum_3clust_lowcond'], colors)
ylabel(['AP count [100ms-1]'])
boxplot_compact({1000*datapool.isi_1clust_lowcond', 1000*datapool.isi_2clust_lowcond', 1000*datapool.isi_3clust_lowcond'}, colors)
ylabel(['ISI [ms]'])
%% Figure S8G
load(fullfile(pwd, 'stat', 'FigureS8G.mat'))
figure
bar(x,mean(a), 'FaceColor',[241,90,36]/255, 'FaceAlpha', 0.6)                
hold on

er = errorbar(x,mean(a),std(a)/sqrt(size(a, 1)));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';                       

%% Figure S9C-D
load(fullfile(pwd, 'stat', 'FigureS9.mat'))
colors = [[0, 0, 0]; [1, 0, 0]];
boxplot_pairwise(p_spike_pool(:,1:2), colors)
boxplot_pairwise(p_burst_pool(:,1:2), colors)
boxplot_pairwise(jitter_pool_strong(:,1:2), colors)
boxplot_pairwise(jitter_pool(:,1:2), colors)
