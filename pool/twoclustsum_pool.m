spike_thres =0.7;
thr = -10;
window_length = 0.1;
dt = 5e-5;
datarange = 19:20;

for i =datarange
    if isempty(data(i).EPSP_singleclust)
        continue
    end
    for n = 1:length(data(i).EPSP_singleclust)
        data(i).amp_single{n} = [];
        data(i).halfwidth_single{n} = [];
        data(i).maxdv_single{n} = [];
        data(i).area_single{n} = [];
        data(i).risetime_single{n} = [];
        if isempty(data(i).EPSP_singleclust{n})
            continue
        end
        for j = 1:size(data(i).EPSP_singleclust{n},2)
            stim_start_new = 500;
            [data(i).amp_single{n}(j),max_amp] = max(abs(data(i).EPSP_singleclust{n}(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_singleclust{n}(stim_start_new, j)));
            max_amp = max_amp+stim_start_new;
            [half,I] = sort(abs((data(i).EPSP_singleclust{n}(:,j)-data(i).EPSP_singleclust{n}(stim_start_new, j))-data(i).amp_single{n}(j)/2));
            b = find(I<max_amp & I>stim_start_new);
            a = find(I>max_amp);
            data(i).halfwidth_single{n}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            data(i).risetime_single{n}(j) = (max_amp-stim_start_new)*dt*1000; %in ms
            data(i).maxdv_single{n}(j) = max(diff(movmean(data(i).EPSP_singleclust{n}(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
            EPSP_mean = data(i).EPSP_singleclust{n}(:,j)-data(i).EPSP_singleclust{n}(stim_start_new, j);
            data(i).area_single{n}(j) = sum(EPSP_mean(501:2501))*dt;
        end
    end
end

for i =datarange
    if isempty(data(i).EPSP_singlespine)
        continue
    end
    for n = 1:length(data(i).EPSP_singlespine)
        data(i).amp_spine{n} = [];
        data(i).halfwidth_spine{n} = [];
        data(i).maxdv_spine{n} = [];
        data(i).area_spine{n} = [];
        data(i).risetime_spine{n} = [];
        if isempty(data(i).EPSP_singlespine{n})
            continue
        end
        for j = 1:size(data(i).EPSP_singlespine{n},2)
            stim_start_new = 500;
            [data(i).amp_spine{n}(j),max_amp] = max(abs(data(i).EPSP_singlespine{n}(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_singlespine{n}(stim_start_new, j)));
            max_amp = max_amp+stim_start_new;
            [half,I] = sort(abs((data(i).EPSP_singlespine{n}(:,j)-data(i).EPSP_singlespine{n}(stim_start_new, j))-data(i).amp_spine{n}(j)/2));
            b = find(I<max_amp & I>stim_start_new);
            a = find(I>max_amp);
            data(i).halfwidth_spine{n}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            data(i).risetime_spine{n}(j) = (max_amp-stim_start_new)*dt*1000; %in ms
            data(i).maxdv_spine{n}(j) = max(diff(movmean(data(i).EPSP_singlespine{n}(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
            EPSP_mean = data(i).EPSP_singlespine{n}(:,j)-data(i).EPSP_singlespine{n}(stim_start_new, j);
            data(i).area_spine{n}(j) = sum(EPSP_mean(501:2501))*dt;
        end
    end
end

for i =datarange
    if isempty(data(i).EPSP_sum)
        continue
    end
    for n = 1:2%length(data(i).EPSP_singleclust)
        data(i).amp_sum{n} = [];
        data(i).halfwidth_sum{n} = [];
        data(i).maxdv_sum{n} = [];
        data(i).area_sum{n} = [];
        data(i).risetime_sum{n} = [];
        data(i).stim_spike_time_lowcond = {};
        if isempty(data(i).EPSP_sum)
            continue
        end
        T = 0:dt:dt*(size(data(i).EPSP_sum,1)-1);
        if n == 1
            for j = 1:7%size(data(i).EPSP_singlespine{n},2)
                stim_start_new = 500;
                [data(i).amp_sum{n}(j),max_amp] = max(abs(data(i).EPSP_sum(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_sum(stim_start_new, j)));
                max_amp = max_amp+stim_start_new;
                [half,I] = sort(abs((data(i).EPSP_sum(:,j)-data(i).EPSP_sum(stim_start_new, j))-data(i).amp_sum{n}(j)/2));
                b = find(I<max_amp & I>stim_start_new);
                a = find(I>max_amp);
                data(i).halfwidth_sum{n}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
                data(i).risetime_sum{n}(j) = (max_amp-stim_start_new)*dt*1000; %in ms
                data(i).maxdv_sum{n}(j) = max(diff(movmean(data(i).EPSP_sum(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
                EPSP_mean = data(i).EPSP_sum(:,j)-data(i).EPSP_sum(stim_start_new, j);
                data(i).area_sum{n}(j) = sum(EPSP_mean(501:2501))*dt;
                dv_over_dt = diff(data(i).EPSP_sum(:, j));
                idx_1 = [];
                for n1 = 2:size(data(i).EPSP_sum,1)
                    if data(i).EPSP_sum(n1, j)> thr && data(i).EPSP_sum(n1-1, j)<= thr
                        idx_1 = [idx_1, n1]; % spike time in s
                    end
                end
                n_spike = length(idx_1);
                [~,c] = findpeaks(dv_over_dt);
                idx = find(dv_over_dt(c)>=spike_thres);
                idx(find(diff(idx)<=2)+1) = [];

                idx_2 = [];
                for n1 = 1:length(idx_1)
                    a = abs(c(idx)-idx_1(n1));
                    I = find(a==min(a));
                    I = I(1);
                    idx_2(n1) = idx(I);
                    idx(I) = [];
                end
                data(i).stim_spike_time_lowcond{j} = T(c(idx_2))-T(stim_start_new);
                data(i).stim_spike_time_lowcond{j} = sort(data(i).stim_spike_time_lowcond{j});
                data(i).stim_spike_time_lowcond{j}(find((data(i).stim_spike_time_lowcond{j}<0) |(data(i).stim_spike_time_lowcond{j}>window_length))) = [];
            end
        else
            for j = 1:8%size(data(i).EPSP_singlespine{n},2)
                stim_start_new = 500;
                [data(i).amp_sum{n}(j),max_amp] = max(abs(data(i).EPSP_sum(stim_start_new:stim_start_new+2000, 7*(n-1)+j)-data(i).EPSP_sum(stim_start_new, 7*(n-1)+j)));
                max_amp = max_amp+stim_start_new;
                [half,I] = sort(abs((data(i).EPSP_sum(:,7*(n-1)+j)-data(i).EPSP_sum(stim_start_new, 7*(n-1)+j))-data(i).amp_sum{n}(j)/2));
                b = find(I<max_amp & I>stim_start_new);
                a = find(I>max_amp);
                data(i).halfwidth_sum{n}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
                data(i).risetime_sum{n}(j) = (max_amp-stim_start_new)*dt*1000; %in ms
                data(i).maxdv_sum{n}(j) = max(diff(movmean(data(i).EPSP_sum(stim_start_new:stim_start_new+1000, 7*(n-1)+j), 10)))/dt/1000; %in V/s
                EPSP_mean = data(i).EPSP_sum(:,7*(n-1)+j)-data(i).EPSP_sum(stim_start_new, 7*(n-1)+j);
                data(i).area_sum{n}(j) = sum(EPSP_mean(501:2501))*dt;
                dv_over_dt = diff(data(i).EPSP_sum(:, j+7*(n-1)));
                idx_1 = [];
                for n1 = 2:size(data(i).EPSP_sum,1)
                    if data(i).EPSP_sum(n1, j+7*(n-1))> thr && data(i).EPSP_sum(n1-1, j+7*(n-1))<= thr
                        idx_1 = [idx_1, n1]; % spike time in s
                    end
                end
                n_spike = length(idx_1);
                [~,c] = findpeaks(dv_over_dt);
                idx = find(dv_over_dt(c)>=spike_thres);
                idx(find(diff(idx)<=2)+1) = [];

                idx_2 = [];
                for n1 = 1:length(idx_1)
                    a = abs(c(idx)-idx_1(n1));
                    I = find(a==min(a));
                    I = I(1);
                    idx_2(n1) = idx(I);
                    idx(I) = [];
                end
                data(i).stim_spike_time_lowcond{j+7*(n-1)} = T(c(idx_2))-T(stim_start_new);
                data(i).stim_spike_time_lowcond{j+7*(n-1)} = sort(data(i).stim_spike_time_lowcond{j+7*(n-1)});
                data(i).stim_spike_time_lowcond{j+7*(n-1)}(find((data(i).stim_spike_time_lowcond{j+7*(n-1)}<0) |(data(i).stim_spike_time_lowcond{j+7*(n-1)}>window_length))) = [];
            end
        end
    end
end

for i = datarange
    if ~isempty(cell2mat(data(i).stim_spike_time_lowcond))
        data(i).onset_lowcond = zeros(1,length(data(i).stim_spike_time_lowcond));
        data(i).APnum_lowcond = zeros(1,length(data(i).stim_spike_time_lowcond));
        data(i).isi_lowcond = cell(1,length(data(i).stim_spike_time_lowcond));
        for j = 1:length(data(i).stim_spike_time_lowcond)            
            if ~isempty(data(i).stim_spike_time_lowcond{j})
                data(i).onset_lowcond(j) = min(data(i).stim_spike_time_lowcond{j});
                data(i).APnum_lowcond(j) = length(data(i).stim_spike_time_lowcond{j});
                data(i).isi_lowcond{j} = diff(data(i).stim_spike_time_lowcond{j});
            end
        end
    end
end
% %%
% dt = 5e-5;
% for i = datarange
%    if (isempty(data(i).EPSP_singleclust))||(data(i).onset == 0)
%        continue
%    end
%    n = 2;
%    stim_start_new = 500;
%    onset_idx = ceil(data(i).onset*1e-3/dt);
%    EPSP_mean = data(i).EPSP_singleclust(:,n)-data(i).EPSP_singleclust(stim_start_new, n);
%    data(i).area_beforeonset = sum(EPSP_mean(stim_start_new+1:stim_start_new+onset_idx))*dt;
%    
% end
%%
datapool.amp = [];
datapool.amp_norm = [];
for i = 1:length(data)
    datapool.amp(i,:) = [data(i).amp_spine{1}(1),data(i).amp_sum{1}, data(i).amp_sum{2}(1:3)];
    datapool.amp_norm(i,:) = datapool.amp(i,:)/max(datapool.amp(i,:));
end
figure
errorbar([1:11],mean(datapool.amp_norm), std(datapool.amp_norm)/sqrt(size(datapool.amp_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'k','MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1] )
ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
datapool.curve_mod = fit([1:11]', mean(datapool.amp_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf]);
hold on
xplot = [0.8:0.1:11.2];
plot(xplot, datapool.curve_mod(xplot), 'k')
xlim([0.8,11.2])
%%
datapool.jitter = [];
datapool.jitter1 = [];
datapool.pburst = [];
datapool.pburst1 = [];
for i = 1:length(data)
    if ~isempty(data(i).jitter)
        a = data(i).jitter([5:7]);
        b = data(i).jitter([9:11]);
        burst  = data(i).p_burst([5:7]);
        burst1  = data(i).p_burst([9:11]);
        datapool.jitter = [datapool.jitter, a(find(a~=0))];
        datapool.jitter1 = [datapool.jitter1, b(find(b~=0))];
        datapool.pburst = [datapool.pburst, burst(find(a~=0))];
        datapool.pburst1 = [datapool.pburst1, burst1(find(b~=0))];
        
    end
end
%%
datapool.onset_lowcond_1clust = [];
datapool.onset_lowcond_2clust = [];
datapool.isi_lowcond_1clust = [];
datapool.isi_lowcond_2clust = [];
datapool.APnum_lowcond_1clust = [];
datapool.APnum_lowcond_2clust = [];

for i = 1:length(data)
    if isempty(cell2mat(data(i).stim_spike_time_lowcond))
        continue
    end
    for j = 8:length(data(i).stim_spike_time_lowcond)   
        if (j<=10)&&(data(i).onset_lowcond(j+5)~=0) &&(data(i).onset_lowcond(j)~=0)
            datapool.onset_lowcond_1clust = [datapool.onset_lowcond_1clust, data(i).onset_lowcond(j)];
            datapool.isi_lowcond_1clust = [datapool.isi_lowcond_1clust, min(data(i).isi_lowcond{j})];
        elseif (j>=13)&&(data(i).onset_lowcond(j)~=0)&&(data(i).onset_lowcond(j-5)~=0)
            datapool.onset_lowcond_2clust = [datapool.onset_lowcond_2clust, data(i).onset_lowcond(j)];
            datapool.isi_lowcond_2clust = [datapool.isi_lowcond_2clust, min(data(i).isi_lowcond{j})];
        end
        if (j<=10)
            datapool.APnum_lowcond_1clust = [datapool.APnum_lowcond_1clust, data(i).APnum_lowcond(j)];
        elseif (j>=13)
            datapool.APnum_lowcond_2clust = [datapool.APnum_lowcond_2clust, data(i).APnum_lowcond(j)];
        end
            
    end
end
%%
datapool.jitter_wostrong_1clust = [];
datapool.jitter_wostrong_2clust = [];
for i = 1:length(data)
    if ~isempty(data(i).jitter)
        if length(data(i).jitter)>7
            if max(cell2mat(data(i).maxdv_single)<2)
                datapool.jitter_wostrong_2clust = [datapool.jitter_wostrong_2clust, data(i).jitter(13:end)];
                datapool.jitter_wostrong_1clust = [datapool.jitter_wostrong_1clust, data(i).jitter(8:10)];
            end
%         else
%             if max(cell2mat(data(i).maxdv_single)<2)
%                 datapool.jitter_wostrong_1clust = [datapool.jitter_wostrong_1clust, data(i).jitter(end-1:end)];
%             end
        end
    end
end
a = find((datapool.jitter_wostrong_1clust==0)|(datapool.jitter_wostrong_2clust==0));
% datapool.jitter_wostrong_1clust(find(datapool.jitter_wostrong_1clust==0)) = [];
% datapool.jitter_wostrong_2clust(find(datapool.jitter_wostrong_2clust==0)) = [];
datapool.jitter_wostrong_1clust(a) = [];
datapool.jitter_wostrong_2clust(a) = [];
%%
count = 0;
datapool.p_burst = [];
datapool.p_spike = [];
datapool.jitter = [];
datapool.onset = [];
for i = 1:length(data)    
    if (isempty(data(i).EPSP_highcond1))||(isempty(data(i).p_burst))
        continue
    end
    count = count+1;
    datapool.p_burst(count,:) = data(i).p_burst;
    datapool.p_spike(count,:) = data(i).p_spike;
    datapool.jitter(count,:) = data(i).jitter;
    data(i).onset = zeros(1,15);
    for j = 1:15
        if isempty(cell2mat(data(i).stim_spike_time(:,j)'))
            continue
        end
        data(i).onset(j) = min(cell2mat(data(i).stim_spike_time(:,j)'))*1e3;
    end
    datapool.onset(count,:) = data(i).onset;
end
%%
count = 0;
datapool.p_burst_1clust = [];
datapool.p_spike_1clust = [];
datapool.jitter_1clust = [];
datapool.onset_1clust = [];
for i = 1:length(data)  
    if isempty(data(i).p_burst)
        continue
    end
    count = count+1;
    datapool.p_burst_1clust(count,:) = data(i).p_burst(1:7);
    datapool.p_spike_1clust(count,:) = data(i).p_spike(1:7);
    datapool.jitter_1clust(count,:) = data(i).jitter(1:7);
    data(i).onset = zeros(1,size(data(i).EPSP_highcond1,2));
    for j = 1:size(data(i).EPSP_highcond1,2)
        if isempty(cell2mat(data(i).stim_spike_time(:,j)'))
            continue
        end
        data(i).onset(j) = min(cell2mat(data(i).stim_spike_time(:,j)'))*1e3;
    end
    datapool.onset_1clust(count,:) = data(i).onset(1:7);
end
%%
for i = 1:8
    data(i).FR1 = 0;
    data(i).FR2 = 0;
    if ~isempty(data(i).EPSP_highcond1)
        data(i).FR1 = length(data(i).spike_time{10})/20;
    end
    if ~isempty(data(i).EPSP_highcond2)
        data(i).FR2 = length(data(i).spike_time2{10})/20;
    end
end
%% 
count = 0;
datapool.amp_nomod = [];
datapool.amp_nomod_1 = [];
datapool.amp_mod = [];
for i = 1:length(data)    
%     if isempty(data(i).EPSP_singleclust{2})
%         continue
%     end
    if isempty(data(i).EPSP_sum)
        continue
    end
    if size(data(i).EPSP_sum, 2)<15
        continue
    end
    count = count + 1;
    if length(data(i).amp_single{2})==7
        if (isempty(data(i).amp_spine))||(isempty(data(i).amp_spine{1})) ||(isempty(data(i).amp_spine{2}))
            datapool.amp_nomod(i,:) = [0, data(i).amp_single{2}];
            
        else
            if length(data(i).amp_spine)==2
                datapool.amp_nomod(i,:) = [data(i).amp_spine{2}(1), data(i).amp_single{2}];
            else
                datapool.amp_nomod(i,:) = [data(i).amp_spine{1}(1), data(i).amp_single{2}];
            end
        end
    else
        datapool.amp_nomod(i,:) = data(i).amp_single{2};
    end
    
    if length(data(i).amp_single{1})==7
        if (isempty(data(i).amp_spine))||(isempty(data(i).amp_spine{1})) ||(isempty(data(i).amp_spine{2}))
            datapool.amp_nomod_1(i,:) = [0, data(i).amp_single{1}];
            
        else
            datapool.amp_nomod_1(i,:) = [data(i).amp_spine{1}(1), data(i).amp_single{1}];
        end
    else
        datapool.amp_nomod_1(i,:) = data(i).amp_single{1};
    end
    
    datapool.amp_mod(i,:) = data(i).amp_sum{2};
    
end
for i = 1:size(datapool.amp_mod, 1)
    datapool.amp_mod_norm(i,:) = datapool.amp_mod(i,:)/max(datapool.amp_mod(i,:));
end
for i = 1:size(datapool.amp_nomod, 1)
    datapool.amp_nomod_norm(i,:) = datapool.amp_nomod(i,:)/max(datapool.amp_nomod(i,:));
    datapool.amp_nomod_norm_1(i,:) = datapool.amp_nomod_1(i,:)/max(datapool.amp_nomod(i,:));
end

%%
datapool.amp_nomod_norm([1,2,4,15, 16],:) = [];
datapool.amp_mod_norm([1,2,4,15, 16],:) = [];
datapool.amp_nomod_norm(8:13,1) = mean(datapool.amp_nomod_norm(1:10,1))*ones(6,1);
datapool.amp_mod_norm([1,2], 4) = mean(datapool.amp_mod_norm(3:13,4))*ones(2,1);

color_nomod = [0,0,0];
color_mod = [1,0,0];
 ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
datapool.curve_mod = fit([1:8]', mean(datapool.amp_mod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude',[4]);
datapool.curve_nomod = fit([1:8]', mean(datapool.amp_nomod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [1]);


data_mod_size = ones(1,8)*size(datapool.amp_mod_norm, 1);
data_mod_size(4) = data_mod_size(4)-2;
data_nomod_size = ones(1,8)*size(datapool.amp_nomod_norm, 1);
data_nomod_size(1) = data_nomod_size(1)-6;

xplot1 = 0:0.1:8.5;
figure
% errorbar([1:8],mean(datapool.amp_mod_norm), std(datapool.amp_mod_norm)/sqrt(size(datapool.amp_mod_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_mod,'MarkerEdgeColor', color_mod,'MarkerFaceColor', [1,1,1] )
errorbar([1:8],mean(datapool.amp_mod_norm), std(datapool.amp_mod_norm)./sqrt(data_mod_size), 'o', 'MarkerSize',10,  'Color', color_mod,'MarkerEdgeColor', color_mod,'MarkerFaceColor', [1,1,1] )
hold on
% errorbar([1:8],mean(datapool.amp_nomod_norm), std(datapool.amp_nomod_norm)/sqrt(size(datapool.amp_nomod_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_nomod,'MarkerEdgeColor', color_nomod,'MarkerFaceColor', [1,1,1] )
errorbar([1:8],mean(datapool.amp_nomod_norm), std(datapool.amp_nomod_norm)./sqrt(data_nomod_size), 'o', 'MarkerSize',10,  'Color', color_nomod,'MarkerEdgeColor', color_nomod,'MarkerFaceColor', [1,1,1] )
hold on
plot(xplot1, datapool.curve_mod(xplot1), 'Color',color_mod)
hold on
plot(xplot1, datapool.curve_nomod(xplot1), 'Color', color_nomod)
xlim([0.8,8.2])
%%
 ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
 
data_nomod_size = ones(1,8)*size(datapool.amp_nomod_norm, 1);
data_nomod_size(1) = data_nomod_size(1)-6;
data_nomod_size_1 = ones(1,8)*size(datapool.amp_nomod_norm_1, 1);
data_nomod_size_1(1) = data_nomod_size(1)-6;
datapool.curve_nomod_1 = fit([1:8]', mean(datapool.amp_nomod_norm_1)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf]);

xplot1 = 0:0.1:8.5;
figure
errorbar([1:8],mean(datapool.amp_nomod_norm), std(datapool.amp_nomod_norm)./sqrt(data_nomod_size), 'o', 'MarkerSize',10,  'Color', [0.5,0.5,0.5],'MarkerEdgeColor', [0.5,0.5,0.5],'MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.amp_nomod_norm_1), std(datapool.amp_nomod_norm_1)./sqrt(data_nomod_size_1), 'o', 'MarkerSize',10,  'Color', 'k','MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1] )
hold on
plot(xplot1, datapool.curve_nomod(xplot1), 'Color',[0.5,0.5,0.5])
hold on
plot(xplot1, datapool.curve_nomod_1(xplot1), 'Color', 'k')
xlim([0.8,8.2])
%%
% N_percolumn =size(datapool.amp_mod, 1)*ones(1,size(datapool.amp_mod, 2));
% for i = 1:size(datapool.amp_mod, 2)
%     idx = find(datapool.amp_mod(:,i)>40);
%     idx1 = find(datapool.amp_mod(:,i)<=40);
%     N_percolumn(i) = N_percolumn(i)-length(idx);
%     for j = 1:length(idx)
%         datapool.amp_mod(idx(j),i) = mean(datapool.amp_mod(idx1,i));
%     end
% end
%
% for i = 1:27
%     data(i).EPSP_singleclust = EPSP_singleclust{i};
%     data(i).amp_single = amp_single_pool{i};
%     data(i).area_single = area_under_curve_pool(i);
%     if i <=length(amp_sum_pool)
%         data(i).amp_multi = amp_sum_pool{i};
%         data(i).EPSP_sum = EPSP_sum{i};
%     end
%     data(i).EPSP_highcond = EPSP_highcond_pool(i,:);
%     data(i).halfwidth_single = halfwidth_single_pool{i};
%     data(i).jitter = jitter_pool(i,:);
%     data(i).maxdv_single = maxdv_single_pool{i};
%     data(i).mean_onset = mean_onset_pool(i,:);
%     data(i).onset = onset_pool(i,:);
%     data(i).p_burst = p_burst_pool(i,:);
%     data(i).p_spike = p_spike_pool(i,:);
%     data(i).risetime_single = risetime_single_pool{i};
%     data(i).spike_time = spike_time_pool(i,:);
%     data(i).stim_spike_time = stim_spike_time_pool{i};
% end
%% samebranch_2dsigmoid
datapool.idx_distmod = [];
datapool.idx_proxmod = [];
for i = [12:85]%1:length(data)
    if data(i).prox_idx==1
        datapool.idx_proxmod = [datapool.idx_proxmod, i];
    else
        datapool.idx_distmod = [datapool.idx_distmod, i];
    end
end

datapool.sum_wodistmod = datapool.amp_nomod(datapool.idx_distmod,:);
datapool.sum_wdistmod = datapool.amp_mod(datapool.idx_distmod,:);
% datapool.sum_wodistmod_norm = datapool.amp_nomod_norm(datapool.idx_distmod,:);
% datapool.sum_wdistmod_norm = datapool.amp_mod_norm(datapool.idx_distmod,:);

datapool.sum_woproxmod = datapool.amp_nomod(datapool.idx_proxmod,:);
datapool.sum_wproxmod = datapool.amp_mod(datapool.idx_proxmod,:);
% datapool.sum_woproxmod_norm = datapool.amp_nomod_norm(datapool.idx_proxmod,:);
% datapool.sum_wproxmod_norm = datapool.amp_mod_norm(datapool.idx_proxmod,:);
for i = 1:size(datapool.sum_wodistmod,1)
    datapool.sum_wodistmod_norm(i,:) = datapool.sum_wodistmod(i,:)/max([datapool.sum_wdistmod(i,:),datapool.sum_wproxmod(i,:)]);
    datapool.sum_woproxmod_norm(i,:) = datapool.sum_woproxmod(i,:)/max([datapool.sum_wdistmod(i,:),datapool.sum_wproxmod(i,:)]);
    datapool.sum_wdistmod_norm(i,:) = datapool.sum_wdistmod(i,:)/max([datapool.sum_wdistmod(i,:),datapool.sum_wproxmod(i,:)]);
    datapool.sum_wproxmod_norm(i,:) = datapool.sum_wproxmod(i,:)/max([datapool.sum_wdistmod(i,:),datapool.sum_wproxmod(i,:)]);
end
datapool.sum_wodistmod_norm([14,26],:) = [];
datapool.sum_woproxmod_norm([14,26],:) = [];
datapool.sum_wdistmod_norm([14,26],:) = [];
datapool.sum_wproxmod_norm([14,26],:) = [];
%%
datapool.area_woproxmod = [];
datapool.area_wproxmod = [];
datapool.area_wodistmod = [];
datapool.area_wdistmod = [];
for i = 1:length(datapool.idx_proxmod)
    datapool.area_wodistmod(i) = data(datapool.idx_distmod(i)).area_sum{1}(end);
    datapool.area_wdistmod(i) = data(datapool.idx_distmod(i)).area_sum{2}(end);
    datapool.area_woproxmod(i) = data(datapool.idx_proxmod(i)).area_sum{1}(end);
    datapool.area_wproxmod(i) = data(datapool.idx_proxmod(i)).area_sum{2}(end);
end
datapool.area_wodistmod([14,26]) = [];
datapool.area_woproxmod([14,26]) = [];
datapool.area_wdistmod([14,26]) = [];
datapool.area_wproxmod([14,26]) = [];
%%
datapool.idx_distmod_Na = [];
datapool.idx_proxmod_Na = [];
for i = 1:length(data)
    if isempty(find((cell2mat(data(i).maxdv_sum)>2)&(cell2mat(data(i).maxdv_sum)<50)))
        continue
    end
    if ~isempty(find(cell2mat(data(i).maxdv_sum)>50))
        continue
    end
    if data(i).prox_idx==1
        datapool.idx_proxmod_Na = [datapool.idx_proxmod_Na, i];
    else
        datapool.idx_distmod_Na = [datapool.idx_distmod_Na, i];
    end
end

datapool.sum_wodistmod_Na = datapool.amp_nomod(datapool.idx_distmod_Na,:);
datapool.sum_wdistmod_Na = datapool.amp_mod(datapool.idx_distmod_Na,:);
datapool.sum_wodistmod_norm_Na = datapool.amp_nomod_norm(datapool.idx_distmod_Na,:);
datapool.sum_wdistmod_norm_Na = datapool.amp_mod_norm(datapool.idx_distmod_Na,:);

datapool.sum_woproxmod_Na = datapool.amp_nomod(datapool.idx_proxmod_Na,:);
datapool.sum_wproxmod_Na = datapool.amp_mod(datapool.idx_proxmod_Na,:);
datapool.sum_woproxmod_norm_Na = datapool.amp_nomod_norm(datapool.idx_proxmod_Na,:);
datapool.sum_wproxmod_norm_Na = datapool.amp_mod_norm(datapool.idx_proxmod_Na,:);

datapool.maxdv_woproxmod = [];
datapool.maxdv_wproxmod = [];
datapool.maxdv_wodistmod = [];
datapool.maxdv_wdistmod = [];
for i = 1:length(datapool.idx_proxmod_Na)
    datapool.maxdv_wodistmod(i) = max(data(datapool.idx_distmod_Na(i)).maxdv_sum{1}(end-2:end));
    datapool.maxdv_wdistmod(i) = max(data(datapool.idx_distmod_Na(i)).maxdv_sum{2}(end-2:end));
    datapool.maxdv_woproxmod(i) = max(data(datapool.idx_proxmod_Na(i)).maxdv_sum{1}(end-2:end));
    datapool.maxdv_wproxmod(i) = max(data(datapool.idx_proxmod_Na(i)).maxdv_sum{2}(end-2:end));
end
%%
dt = 5e-5;
datapool.amp_linearsum = [];
datapool.area_linearsum = [];
datapool.maxdv_linearsum = [];
i_range = 12:2:85;
for i = 1:length(i_range)
    linearsum = data(i_range(i)).EPSP_singleclust{1}(:,end) + data(i_range(i)).EPSP_singleclust{2}(:,end);
    stim_start_new = 500;
    [datapool.amp_linearsum(i),max_amp] = max(abs(linearsum(stim_start_new:stim_start_new+2000)-linearsum(stim_start_new)));
    max_amp = max_amp+stim_start_new;
    [half,I] = sort(abs((linearsum-linearsum(stim_start_new))-datapool.amp_linearsum(i)/2));
    b = find(I<max_amp & I>stim_start_new);
    a = find(I>max_amp);
    EPSP_mean = linearsum-linearsum(stim_start_new);
    datapool.area_linearsum(i) = sum(EPSP_mean(501:2501))*dt;
    if (~isempty(find(datapool.idx_proxmod_Na==i_range(i))))||(~isempty(find(datapool.idx_proxmod_Na==i_range(i)+1)))
        datapool.maxdv_linearsum = [datapool.maxdv_linearsum, max(diff(movmean(linearsum(stim_start_new:stim_start_new+1000, :), 10)))/dt/1000]; %in V/s
    end
end
datapool.area_linearsum([14,26]) = [];
datapool.amp_linearsum([14,26]) = [];

%%
figure
errorbar([1:8],mean(datapool.sum_wdistmod_norm), std(datapool.sum_wdistmod_norm)/sqrt(size(datapool.sum_wdistmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'r','MarkerEdgeColor', 'r','MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.sum_wodistmod_norm), std(datapool.sum_wodistmod_norm)/sqrt(size(datapool.sum_wodistmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'k','MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1] )

figure
errorbar([1:8],mean(datapool.sum_wproxmod_norm), std(datapool.sum_wproxmod_norm)/sqrt(size(datapool.sum_wproxmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'b','MarkerEdgeColor', 'b','MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.sum_woproxmod_norm), std(datapool.sum_woproxmod_norm)/sqrt(size(datapool.sum_woproxmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'k','MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1] )

%%
color_prox = [147,39,143]/255;
color_dist = [0,164,79]/255;
 ft = fittype('d+b/(1+exp(-(x-a)/k))', 'independent', 'x');
datapool.curve_wdistmod = fit([1:8]', mean(datapool.sum_wdistmod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [7]);
datapool.curve_wodistmod = fit([1:8]', mean(datapool.sum_wodistmod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [1,2]);

xplot1 = 1:0.1:8.5;
figure
errorbar([1:8],mean(datapool.sum_wdistmod_norm), std(datapool.sum_wdistmod_norm)/sqrt(size(datapool.sum_wdistmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'r','MarkerEdgeColor', 'r','MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.sum_wodistmod_norm), std(datapool.sum_wodistmod_norm)/sqrt(size(datapool.sum_wodistmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_prox,'MarkerEdgeColor', color_prox,'MarkerFaceColor', [1,1,1] )
hold on
plot(xplot1, datapool.curve_wodistmod(xplot1), 'Color',color_prox)
hold on
plot(xplot1, datapool.curve_wdistmod(xplot1), 'r')
xlim([1.8,8.2])

%%
datapool.curve_wproxmod = fit([1:8]', mean(datapool.sum_wproxmod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf]);
datapool.curve_woproxmod = fit([1:8]', mean(datapool.sum_woproxmod_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [1]);

figure
errorbar([1:8],mean(datapool.sum_wproxmod_norm), std(datapool.sum_wproxmod_norm)/sqrt(size(datapool.sum_wproxmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', 'b','MarkerEdgeColor', 'b','MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.sum_woproxmod_norm), std(datapool.sum_woproxmod_norm)/sqrt(size(datapool.sum_woproxmod_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_dist,'MarkerEdgeColor', color_dist,'MarkerFaceColor', [1,1,1] )
hold on
plot(xplot1, datapool.curve_woproxmod(xplot1), 'Color', color_dist)
hold on
plot(xplot1, datapool.curve_wproxmod(xplot1), 'b')
xlim([1.8,8.2])
%%
for i = 1:size(datapool.sum_wproxmod_norm,1)
    curve_wproxmod = fit([1:8]', datapool.sum_wproxmod_norm(i,:)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, 2], 'Lower', [0, -inf, 0, 0]);
    curve_wdistmod = fit([1:8]', datapool.sum_wdistmod_norm(i,:)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, 2], 'Lower', [0, -inf, 0, 0]);
    curve_woproxmod = fit([1:8]', datapool.sum_woproxmod_norm(i,:)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, 2], 'Lower', [0, -inf, 0, 0], 'Exclude',[1]);
    curve_wodistmod = fit([1:8]', datapool.sum_wodistmod_norm(i,:)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, 2], 'Lower', [0, -inf, 0, 0], 'Exclude', [1]);
    datapool.sigmoid_k(i,:) = [curve_wproxmod.k,curve_woproxmod.k, curve_wdistmod.k, curve_wodistmod.k];
   
end
datapool.sigmoid_k_exclude = datapool.sigmoid_k;
idx_exclude = [];
for i = 1:size(datapool.sigmoid_k,1)
    if ~isempty(find(abs(datapool.sigmoid_k(i,:)-2)<=0.1, 1))
        idx_exclude = [idx_exclude,i];
    end
end
datapool.sigmoid_k_exclude(idx_exclude,:) = [];


figure
boxplot(datapool.sigmoid_k_exclude,'Color', [[0,0,1];color_dist;[1,0,0];color_prox],'PlotStyle','compact')
%%
% count = 0;
% datapool.sum_prox = [];
% datapool.sum_dist = [];
% for i = 1:length(data)
%     if isempty(data(i).EPSP_singleclust)
%         continue
%     end
%     if isempty(data(i).EPSP_singleclust{1})
%         continue
%     end
%     if isempty(data(i).EPSP_singleclust{2})
%         continue
%     end
%     count = count +1;
%     datapool.sum_prox(count,:) = data(i).amp_single{data(i).prox_idx}(:,end-6:end);
%     datapool.sum_dist(count,:) = data(i).amp_single{data(i).dist_idx}(:,end-6:end);
% end
% for i = 1:size(datapool.sum_prox, 1)
%     datapool.sum_prox_norm(i,:) = datapool.sum_prox(i,:)/max(datapool.sum_prox(i,:));
%     datapool.sum_dist_norm(i,:) = datapool.sum_dist(i,:)/max(datapool.sum_dist(i,:));
% end
datapool.sum_prox = datapool.sum_wodistmod;
datapool.sum_dist = datapool.sum_woproxmod;
datapool.sum_prox_norm = [];
datapool.sum_dist_norm = [];
for i = 1:size(datapool.sum_prox, 1)
    datapool.sum_prox_norm(i,:) = datapool.sum_prox(i,:)/max(datapool.sum_prox(i,:));
    
end
for i = 1:size(datapool.sum_dist, 1)
    datapool.sum_dist_norm(i,:) = datapool.sum_dist(i,:)/max(datapool.sum_dist(i,:));
end

figure
errorbar([1:8],mean(datapool.sum_prox_norm), std(datapool.sum_prox_norm)/sqrt(size(datapool.sum_prox_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_prox,'MarkerEdgeColor', color_prox,'MarkerFaceColor', [1,1,1] )
hold on
errorbar([1:8],mean(datapool.sum_dist_norm), std(datapool.sum_dist_norm)/sqrt(size(datapool.sum_dist_norm, 1)), 'o', 'MarkerSize',10,  'Color', color_dist,'MarkerEdgeColor', color_dist,'MarkerFaceColor', [1,1,1] )

datapool.curve_prox = fit([1:8]', mean(datapool.sum_prox_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [1,2]);
datapool.curve_dist = fit([1:8]', mean(datapool.sum_dist_norm)',ft, 'Start', [4, 1, 0,1], 'Upper', [8, inf, 1, inf], 'Lower', [0, -inf, 0, -inf], 'Exclude', [1]);
xplot1 = 1:0.1:8.5;
hold on
plot(xplot1, datapool.curve_prox(xplot1), 'Color', color_prox, 'Linewidth', 1)
plot(xplot1, datapool.curve_dist(xplot1), 'Color', color_dist, 'Linewidth', 1)
xlim([1.8,8.2])

%% discrimination between clusters
% datarange = 1:length(data);
% biuld kernel
win_length = 0.1; %100ms
bin_size = 1e-4; %5ms
edges = 0:bin_size:win_length;
dt1 = bin_size*1e3;
epsp_kernel = zeros(1,150);
A = 0.005;
B = 0.005;
tau1 = 0.1;
tau2 = 2;
tp = (tau1*tau2)/(tau2-tau1)*log(tau2/tau1);
factor = 1/(-exp(-tp/tau1)+exp(-tp/tau2));
for i = 1:length(epsp_kernel)-1
    A = A+dt1*(-A/tau1);
    B = B+dt1*(-B/tau2);
    epsp_kernel(i+1) = B-A;%epsp_kernel(i)+dt*(-epsp_kernel(i)/tau2 + xepsp*(1-epsp_kernel(i)));
end
w = epsp_kernel/sum(epsp_kernel);
for i = datarange
    for j = 1:5
        data(i).stim_spike_time_matrix{j} = zeros(size(data(i).stim_spike_time, 2), length(edges)-1);
        data(i).stim_spike_time_kernel{j} = zeros(size(data(i).stim_spike_time, 2), length(edges)-1);
    end
    for j = 1:size(data(i).stim_spike_time, 1)
        for n = 1:size(data(i).stim_spike_time, 2)
            if isempty(data(i).stim_spike_time{j,n})
                continue
            else
                data(i).stim_spike_time_matrix{j}(n,:) = histcounts(data(i).stim_spike_time{j,n}, edges);
                a = conv(data(i).stim_spike_time_matrix{j}(n,:),w);
                data(i).stim_spike_time_kernel{j}(n,:) = a(1:size(data(i).stim_spike_time_matrix{j}(n,:), 2));
            end
        end
    end            
end
%% 
% decode, build stim and spontaneous kernels

win_length = 0.1; %100ms
bin_size = 1e-4; %5ms
edges = 0:bin_size:win_length;
dt1 = bin_size*1e3;
epsp_kernel = zeros(1,150);
A = 0.005;
B = 0.005;
tau1 = 0.1;
tau2 = 2;
tp = (tau1*tau2)/(tau2-tau1)*log(tau2/tau1);
factor = 1/(-exp(-tp/tau1)+exp(-tp/tau2));
for i = 1:length(epsp_kernel)-1
    A = A+dt1*(-A/tau1);
    B = B+dt1*(-B/tau2);
    epsp_kernel(i+1) = B-A;%epsp_kernel(i)+dt*(-epsp_kernel(i)/tau2 + xepsp*(1-epsp_kernel(i)));
end
w = epsp_kernel/sum(epsp_kernel);


for i = datarange
    data(i).strong_idx = [];
    for j = 1:length(data(i).maxdv_single)
    if ~isempty(find((data(i).maxdv_single{j}>2)&(data(i).maxdv_single{j}<50)))
        data(i).strong_idx = [data(i).strong_idx,j];
    end
    end
    
    if isempty(data(i).stim_spike_time)
        continue
    end
    for j = 1:size(data(i).stim_spike_time, 2)
        data(i).stim_spike_time_matrix{j} = zeros(size(data(i).stim_spike_time, 1), length(edges)-1);
        data(i).stim_spike_time_kernel{j} = zeros(size(data(i).stim_spike_time, 1), length(edges)-1);
    end

    for j = 1:size(data(i).stim_spike_time, 2)
        for n = 1:size(data(i).stim_spike_time, 1)
            if isempty(cell2mat(data(i).stim_spike_time(:,j)'))
                continue
            else
                data(i).stim_spike_time_matrix{j}(n,:) = histcounts(data(i).stim_spike_time{n,j}, edges);
                a = conv(data(i).stim_spike_time_matrix{j}(n,:),w);
                data(i).stim_spike_time_kernel{j}(n,:) = a(1:size(data(i).stim_spike_time_matrix{j}(n,:), 2));
            end
        end
    end 

    data(i).spon_spike_time = [];
    trial_break = 0;
    for j = 1:length(data(i).spike_time)
        spike_time = data(i).spike_time{j};
        orig_spike_time = data(i).spike_time{j};
        flag = 0;
        stimtrig_idx = [];
        if j >length(data(i).stim_start_time)
            stim_trig_idx = [];
        else
            for m = 1:length(data(i).stim_start_time{j})
                stimtrig_idx = [stimtrig_idx, find((orig_spike_time-data(i).stim_start_time{j}(m)<=win_length)&(orig_spike_time-data(i).stim_start_time{j}(m)>0))];
                if ~isempty(find((orig_spike_time-data(i).stim_start_time{j}(m)<=win_length)&(orig_spike_time-data(i).stim_start_time{j}(m)>0)))
                    spike_time(stimtrig_idx(end)+1:end) = spike_time(stimtrig_idx(end)+1:end)-win_length;
                    flag = flag+1;
                end
            end
        end
        spike_time(stimtrig_idx) = [];
        data(i).spon_spike_time = [data(i).spon_spike_time, spike_time + trial_break];
        trial_break = trial_break + 10-win_length*flag;
    end

    if isempty(data(i).spon_spike_time)
        data(i).spon_spike_time_matrix = zeros(1,1e4);
        data(i).spon_spike_time_kernel = zeros(1,1e4);
        continue
    end
    spon_edges = 0:bin_size:ceil(data(i).spon_spike_time(end)/win_length)*win_length;
    data(i).spon_spike_time_matrix = histcounts(data(i).spon_spike_time, spon_edges);
    if length(data(i).spon_spike_time_matrix) < 1e4
        data(i).spon_spike_time_matrix = [data(i).spon_spike_time_matrix, zeros(1, 1e4-length(data(i).spon_spike_time_matrix))];
    end
    a = conv(data(i).spon_spike_time_matrix,w);
    data(i).spon_spike_time_kernel = a(1:length(data(i).spon_spike_time_matrix));
        
end

for i = datarange
    data(i).spon_FR = [];
    if isempty(data(i).stim_spike_time)
        continue
    end
    data(i).spon_FR = length(find(data(i).spon_spike_time_matrix~=0))/(length(data(i).spon_spike_time_matrix)*bin_size);
end

%% discrimination between clusters with diff.driving
datapool.n_trials = 10;
for i = datarange
    data(i).linear_regression_R = [];
    data(i).linear_regression_error = [];
    data(i).linear_regression_R_shuffle_clust = [];
    data(i).linear_regression_error_shuffle_clust = [];
end
for i = 12:2:85%datarange
    if i<12
        continue
    end
    if i == 62
        continue
    end
    if isempty(data(i).stim_spike_time)
        continue
    end
    
    if data(i).prox_idx == 1
        A = [cell2mat(data(i+1).stim_spike_time_kernel(7:end)');cell2mat(data(i).stim_spike_time_kernel(7:end)')];
    else
        A = [cell2mat(data(i).stim_spike_time_kernel(7:end)');cell2mat(data(i+1).stim_spike_time_kernel(7:end)')];
    end
    [data(i).linear_regression_R, data(i).linear_regression_error] = multicluster_linear_regression(A, datapool.n_trials);
    [data(i).linear_regression_R_shuffle_clust, data(i).linear_regression_error_shuffle_clust] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
    for n = 1:9
        [R_shuffle_clust, error_shuffle_clust] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
        data(i).linear_regression_R_shuffle_clust = data(i).linear_regression_R_shuffle_clust + R_shuffle_clust;
        data(i).linear_regression_error_shuffle_clust = data(i).linear_regression_error_shuffle_clust + error_shuffle_clust;
    end
    data(i).linear_regression_R_shuffle_clust = data(i).linear_regression_R_shuffle_clust/10;
    data(i).linear_regression_error_shuffle_clust = data(i).linear_regression_error_shuffle_clust/10;
    disp(sprintf('cell number %d done', i));
end

%% pool of decode between combination of partial cluster

datapool.R_pool = [];
datapool.R_shuffle_pool = [];
count = 1;
for i = 1:length(data)
    if isempty(data(i).linear_regression_R)
        continue
    end
    datapool.R_pool(count,:) = diag(data(i).linear_regression_R,9)';
    datapool.R_shuffle_pool(count,:) = diag(data(i).linear_regression_R_shuffle_clust,9)';
    count = count+1;
end
