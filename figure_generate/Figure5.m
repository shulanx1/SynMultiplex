%% Figure 5A
dt = 5e-5;
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
idx_cell = 12;  
idx_trials = [1, 2];  
idx_cond = 3;

t = -500*dt:dt:dt*(size(data(data_decode(idx_cell).idx(idx_trials(1))).EPSP_highcond{idx_cond}, 1)-1-500);
for i = 1:length(idx_trials)
    figure
    subplot(131)
    for j = 1:size(data(data_decode(idx_cell).idx(idx_trials(i))).EPSP_highcond{idx_cond}, 2)
        plot(t, data(data_decode(idx_cell).idx(idx_trials(i))).EPSP_highcond{idx_cond}(:,j)-50*(j-1), 'k')
        hold on
    end
    xlim([-0.01,0.1])
    subplot(132)
    raster_plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time(idx_cond,:))
    subplot(133)
    for j = 1:size(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time_kernel{idx_cond}, 1)
        plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time_kernel{idx_cond}(j,:)-0.05*(j-1), 'k')
        hold on
    end
end
%% Figure 5B
addpath(genpath(fullfile(pwd, 'GeneSetAnalysisMatlab')))
figure, imagesc(data_decode(idx_cell).correlation_matrix{idx_cond}(10*(idx_trials(1)-1)+1:10*(idx_trials(end)),10*(idx_trials(1)-1)+1:10*(idx_trials(end))))
cmap = custom_cmap('redbluedark');
colormap(cmap), colorbar()
%% Figure 5D
A = [data(data_decode(idx_cell).idx(idx_trials(1))).stim_spike_time_kernel{idx_cond};data(data_decode(idx_cell).idx(idx_trials(2))).stim_spike_time_kernel{idx_cond}];
[R, error] = multicluster_linear_regression(A, 10, 1);
%% Figure 5E
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode_ap5.mat'))
figure
colors = [[247,147,30];[128,128,128];[0, 145, 69];[128,128,128]]/255;
boxplot_pairwise({cell2mat(datapool.R_betweenclust_pool(2:5)), cell2mat(datapool_ap5.R_betweenclust_pool(2:5))}, colors)
%% Figure 5F
figure
colors = [[247,147,30];[128,128,128]]/255;
boxplot_pairwise({datapool.R_betweenclust_pool{1}, datapool.R_betweenclust_pool{2}, datapool.R_betweenclust_pool{4} }, colors)
boxplot_pairwise({cell2mat(datapool.R_betweenclust_pool(2:3)), cell2mat(datapool.R_betweenclust_pool(4:5))}, colors)
%% Figure 5H
x = [1:2:9];
colors = [[247,147,30];[128,128,128]]/255;
for i = 1:5
    for j = 1:2
        y(j,i) = mean(datapool.R_FR{i}(:,j));
        ymse(j,i) = std(datapool.R_FR{i}(:,j))/sqrt(size(datapool.R_FR{i}, 1));
    end
end
figure
for i = 1:2
    color_idx = mod(i, size(colors, 1));
    if color_idx == 0
        color_idx = size(colors, 1);
    end
    plot(x, y(i,:), 'Color', colors(color_idx,:),'Linewidth', 2)
    hold on
    plot(x, y(i,:) + ymse(i,:), 'Color', colors(color_idx,:),'Linewidth', 0.5)
    plot(x, y(i,:) - ymse(i,:), 'Color', colors(color_idx,:),'Linewidth', 0.5)
end
%% FigureS10D
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode_ap5.mat'))
figure
colors = [[247,147,30];[128,128,128];[0,0,0];[128,128,128];[0, 145, 69];[128,128,128];[0,0,0];[128,128,128]]/255;
boxplot_pairwise({cell2mat(datapool.R_clustandspon_pool(2:5)), cell2mat(datapool_ap5.R_clustandspon_pool(2:5))}, colors)