%% Figure 6A-C
dt = 5e-5;
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
idx_cell = 19;  
idx_trials = [12,13,14,16];  
idx_cond = 2;
t = -500*dt:dt:dt*(size(data(data_decode(idx_cell).idx(idx_trials(1))).EPSP_singleclust, 1)-1-500);
figure
plot(t, data(data_decode(idx_cell).idx(idx_trials(1))).EPSP_singleclust)
figure
for i = 1:length(idx_trials)
    subplot(1,length(idx_trials), i)
    raster_plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time(idx_cond,:))
end
figure,imagesc(data_decode(idx_cell).linear_regression_R{idx_cond}(idx_trials, idx_trials))
cmap = custom_cmap('redbluedark');
colormap(cmap), colorbar(), caxis([0,0.5])
%% Figure 6D
datapool.R_strong = cell(1,5);
datapool.R_weak = cell(1,5);
datapool.R_strongweak = cell(1,5);
for i = 1:length(datapool.samecell_idx)
    for j = 1:5
        for n = 1:length(data_decode(i).clust)-1
            clust = data_decode(i).clust{n};
            clust(find(clust==999)) = [];
            if length(clust)>=4
                continue;
            end
            if size(data_decode(i).linear_regression_R{j}, 1)<n
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)<data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(1,2)
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)==0
                continue
            end
            for m = n+1:length(data_decode(i).clust)
                clust2 = data_decode(i).clust{n};
                clust2(find(clust2==999)) = [];
                if length(clust2)~=length(clust)
                    continue;
                end
                if size(data_decode(i).linear_regression_R{j}, 1)<max(n,m)
                    continue
                end
                if ((data_decode(i).strong_num(n)>=1)&&(data_decode(i).strong_num(m)<1))||((data_decode(i).strong_num(n)<1)&&(data_decode(i).strong_num(m)>=1))
                    datapool.R_strongweak{j} = [datapool.R_strongweak{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((data_decode(i).strong_num(n)>=1)&&(data_decode(i).strong_num(m)>=1))
                     datapool.R_strong{j} = [datapool.R_strong{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((data_decode(i).strong_num(n)<1)&&(data_decode(i).strong_num(m)<1))
                    datapool.R_weak{j} = [datapool.R_weak{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                end
            end
        end

    end
end
figure
colors = [[165,91,90];[128,128,128];[3,87,134];[128,128,128]]/255;
boxplot_pairwise({cell2mat(datapool.R_strong(2:5)'), cell2mat(datapool.R_strongweak(2:5)')}, colors)
%% Figure 6E-H
dt = 5e-5;
load(fullfile(pwd, 'datasets', 'highcond_multibranch_decode.mat'))
idx_cell = 20;  
idx_trials = [1:8];  
idx_cond = 2;
t = -500*dt:dt:dt*(size(data(data_decode(idx_cell).idx(idx_trials(1))).EPSP_singleclust, 1)-1-500);
figure
for i = 1:length(idx_trials)
    subplot(1,length(idx_trials), i)
    raster_plot(data(data_decode(idx_cell).idx(idx_trials(i))).stim_spike_time(idx_cond,:))
end
figure,imagesc(data_decode(idx_cell).linear_regression_R{idx_cond}(idx_trials, idx_trials))
cmap = custom_cmap('redbluedark');
colormap(cmap), colorbar(),caxis([0,0.8])
%% Figure 6I
figure
colors = [[165,91,90];[128,128,128];[3,87,134];[128,128,128]]/255;
boxplot_pairwise({cell2mat(datapool.R_strong_4clust(2:5)'), cell2mat(datapool.R_weak_4clust(2:5)')}, colors)