% compute width etc. for single cluster EPSP
datarange = 64:69;%1:length(data);

%%
for i = datarange
    if isempty(data(i).EPSP_singleclust)
        continue
    end
    data(i).amp_single = [];
    data(i).halfwidth_single = [];
    data(i).maxdv_single = [];
    data(i).area_single = [];
    data(i).risetime_single = [];
    if isempty(data(i).EPSP_singleclust)
        continue
    end
    for j = 1:size(data(i).EPSP_singleclust,2)
        stim_start_new = 500;
        [data(i).amp_single(j),max_amp] = max(abs(data(i).EPSP_singleclust(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_singleclust(stim_start_new, j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs((data(i).EPSP_singleclust(:,j)-data(i).EPSP_singleclust(stim_start_new, j))-data(i).amp_single(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        data(i).halfwidth_single(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
        data(i).risetime_single(j) = (max_amp-stim_start_new)*dt*1000; %in ms
        data(i).maxdv_single(j) = max(diff(movmean(data(i).EPSP_singleclust(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
        EPSP_mean = data(i).EPSP_singleclust(:,j)-data(i).EPSP_singleclust(stim_start_new, j);
        data(i).area_single(j) = sum(EPSP_mean(501:2501))*dt;
    end

end

% compute width etc. and spike time for sum EPSP
spike_thres =0.5;
thr = -10;
window_length = 0.1;
dt = 5e-5;
for i = datarange
    if isempty(data(i).EPSP_sum)
        continue
    end
        data(i).amp_sum = [];
        data(i).halfwidth_sum = [];
        data(i).maxdv_sum = [];
        data(i).area_sum = [];
        data(i).risetime_sum = [];
        data(i).stim_spike_time_lowcond = {};
        if isempty(data(i).EPSP_sum)
            continue
        end
        T = 0:dt:dt*(size(data(i).EPSP_sum,1)-1);
        for j = 1:size(data(i).EPSP_sum,2)
            stim_start_new = 500;
            [data(i).amp_sum(j),max_amp] = max(abs(data(i).EPSP_sum(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_sum(stim_start_new, j)));
            max_amp = max_amp+stim_start_new;
            [half,I] = sort(abs((data(i).EPSP_sum(:,j)-data(i).EPSP_sum(stim_start_new, j))-data(i).amp_sum(j)/2));
            b = find(I<max_amp & I>stim_start_new);
            a = find(I>max_amp);
            data(i).halfwidth_sum(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            data(i).risetime_sum(j) = (max_amp-stim_start_new)*dt*1000; %in ms
            data(i).maxdv_sum(j) = max(diff(movmean(data(i).EPSP_sum(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
            EPSP_mean = data(i).EPSP_sum(:,j)-data(i).EPSP_sum(stim_start_new, j);
            data(i).area_sum(j) = sum(EPSP_mean(501:2501))*dt;
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
 
end
%% compute onset, AP number, isi for high cond state
for i = datarange
    data(i).mean_onset = [];
    data(i).median_onset = [];
    data(i).onset = [];
    data(i).AP_num = [];
    data(i).isi = {};
    if isempty(data(i).stim_spike_time)
        continue
    end
    data(i).mean_onset = zeros(1,5);
    data(i).onset = zeros(1,5);
    data(i).median_onset = zeros(1,5);
    data(i).AP_num = zeros(1,5);
    for j = 1:size(data(i).stim_spike_time,1)
        data(i).isi{j} = [];
        if isempty(cell2mat(data(i).stim_spike_time(j,:)))
            continue
        else
            a = [];
            for n = 1:size(data(i).stim_spike_time,2)
                if ~isempty(data(i).stim_spike_time{j,n})
                    a = [a,min(data(i).stim_spike_time{j,n})];
                end
            end
            data(i).mean_onset(j) = mean(a)*1e3;
            data(i).median_onset(j) = median(a)*1e3;
            data(i).onset(j) = min(cell2mat(data(i).stim_spike_time(j,:)))*1e3;
            data(i).AP_num(j) = length(cell2mat(data(i).stim_spike_time(j,:)))/10;
            for n = 1:size(data(i).stim_spike_time,2)
                if length(data(i).stim_spike_time{j,n})<=1
                    continue
                else
                    data(i).isi{j} = [data(i).isi{j}, diff(data(i).stim_spike_time{j,n})];
                end
            end
        end
    end
end
for i = datarange
    
    if isempty(data(i).stim_spike_time)
        continue
    end
    data(i).p_burst_true = zeros(1,5);
    for j = 1:size(data(i).stim_spike_time,1)
        if isempty(cell2mat(data(i).stim_spike_time(j,:)))
            continue
        else
            p_burst = 0;
            for n = 1:size(data(i).stim_spike_time,2)
                if length(data(i).stim_spike_time{j,n})<=1
                    continue
                else                    
                    if min(diff(data(i).stim_spike_time{j,n})*1e3)<= 20
                        p_burst = p_burst + 1;
                    end
                end
            end
            data(i).p_burst_true(j) = p_burst/size(data(i).stim_spike_time,2);
        end
    end
end
%% pool p burst etc.
count = 0;
datapool.p_burst = [];
datapool.p_spike = [];
datapool.jitter = [];
datapool.onset = [];
datapool.mean_onset = [];
datapool.median_onset = [];
datapool.isi_lowcond = {};
datapool.isi_highcond = {};
datapool.AP_num = [];
for i = 1:length(data)   
    if isempty(data(i).EPSP_highcond)
        continue
    end
    count = count+1;
    datapool.AP_num(count, 1:length(data(i).AP_num)) = data(i).AP_num;
    datapool.p_burst(count,1:length(data(i).p_burst)) = data(i).p_burst;
    datapool.p_burst_true(count,1:length(data(i).p_burst_true)) = data(i).p_burst_true;
    datapool.p_spike(count,1:length(data(i).p_spike)) = data(i).p_spike;
    datapool.jitter(count,1:length(data(i).jitter)) = data(i).jitter;
    datapool.onset(count,1:length(data(i).onset)) = data(i).onset;
    datapool.mean_onset(count,1:length(data(i).mean_onset)) = data(i).mean_onset;
    datapool.median_onset(count,1:length(data(i).median_onset)) = data(i).median_onset;
    datapool.isi_lowcond(count) = data(i).isi(1);
    if length(data(i).isi)>1
        datapool.isi_highcond(count) = data(i).isi(2);
    end
end
%% pool together for all conditions

datapool.idx_strong = [];
datapool.idx_weak = [];
datapool.idx_burst_lowcond = [];
datapool.idx_single_lowcond = [];
datapool.idx_burst_highcond = [];
datapool.idx_single_highcond = [];
datapool.idx_burst_highcond_highI = [];
datapool.idx_single_highcond_highI = [];
for i = 1:length(data)
    clust = data(i).clust;
    clust(find(clust==999)) = [];
    if length(clust)<2
        continue
    end
    if ~isempty(find((data(i).maxdv_single>2)&(data(i).maxdv_single<50)))
        datapool.idx_strong = [datapool.idx_strong,i];
    else
        datapool.idx_weak = [datapool.idx_weak,i];
    end
    if (data(i).p_burst(1)>0.3)&(data(i).jitter(1)~=0)
        datapool.idx_burst_lowcond = [datapool.idx_burst_lowcond, i];
    elseif (data(i).p_burst(1)<=0.3)&&(data(i).jitter(1)~=0)
        datapool.idx_single_lowcond = [datapool.idx_single_lowcond, i];
    end
    if (data(i).p_burst_true(2)>0.3)&&(data(i).jitter(2)~=0)
        datapool.idx_burst_highcond = [datapool.idx_burst_highcond, i];
    elseif (data(i).p_burst_true(2)<=0.3)&&(data(i).jitter(2)~=0)
        datapool.idx_single_highcond = [datapool.idx_single_highcond, i];
    end
    
end
for i = 1:length(data)
    clust = data(i).clust;
    clust(find(clust==999)) = [];
    if length(clust)<2
        continue
    end
    if length(data(i).p_burst)<4
        continue
    end
    if length(data(i).jitter)<4
        continue
    end
    if (data(i).p_burst_true(4)>0.3)&&(data(i).jitter(4)~=0)
        datapool.idx_burst_highcond_highI = [datapool.idx_burst_highcond_highI, i];
    elseif (data(i).p_burst_true(4)<=0.3)&&(data(i).jitter(4)~=0)
        datapool.idx_single_highcond_highI = [datapool.idx_single_highcond_highI, i];
    end
end

datapool.jitter_strong_lowcond = datapool.jitter(datapool.idx_strong,1);
datapool.jitter_strong_lowcond(find(datapool.jitter_strong_lowcond==0))= [];
datapool.jitter_weak_lowcond = datapool.jitter(datapool.idx_weak,1);
datapool.jitter_weak_lowcond(find(datapool.jitter_weak_lowcond==0))= [];
datapool.jitter_strong_highcond = datapool.jitter(datapool.idx_strong,2);
datapool.jitter_strong_highcond(find(datapool.jitter_strong_highcond==0))= [];
datapool.jitter_strong_highcond_highI = datapool.jitter(datapool.idx_strong,4);
datapool.jitter_strong_highcond_highI(find(datapool.jitter_strong_highcond_highI==0))= [];
datapool.jitter_weak_highcond = datapool.jitter(datapool.idx_weak,2);
datapool.jitter_weak_highcond(find(datapool.jitter_weak_highcond==0))= [];
datapool.jitter_weak_highcond_highI = datapool.jitter(datapool.idx_weak,4);
datapool.jitter_weak_highcond_highI(find(datapool.jitter_weak_highcond_highI ==0))= [];

datapool.mean_onset_strong_lowcond = datapool.mean_onset(datapool.idx_strong,1);
datapool.mean_onset_strong_lowcond(find(datapool.mean_onset_strong_lowcond==0))= [];
datapool.mean_onset_weak_lowcond = datapool.mean_onset(datapool.idx_weak,1);
datapool.mean_onset_weak_lowcond(find(datapool.mean_onset_weak_lowcond==0))= [];
datapool.mean_onset_strong_highcond = datapool.mean_onset(datapool.idx_strong,2);
datapool.mean_onset_strong_highcond(find(datapool.mean_onset_strong_highcond==0))= [];
datapool.mean_onset_weak_highcond = datapool.mean_onset(datapool.idx_weak,2);
datapool.mean_onset_weak_highcond(find(datapool.mean_onset_weak_highcond==0))= [];
datapool.mean_onset_strong_highcond_highI = datapool.mean_onset(datapool.idx_strong,4);
datapool.mean_onset_strong_highcond_highI(find(datapool.mean_onset_strong_highcond_highI==0))= [];
datapool.mean_onset_weak_highcond_highI = datapool.mean_onset(datapool.idx_weak,4);
datapool.mean_onset_weak_highcond_highI(find(datapool.mean_onset_weak_highcond_highI==0))= [];

datapool.median_onset_strong_lowcond = datapool.median_onset(datapool.idx_strong,1);
datapool.median_onset_strong_lowcond(find(datapool.median_onset_strong_lowcond==0))= [];
datapool.median_onset_weak_lowcond = datapool.median_onset(datapool.idx_weak,1);
datapool.median_onset_weak_lowcond(find(datapool.median_onset_weak_lowcond==0))= [];
datapool.median_onset_strong_highcond = datapool.median_onset(datapool.idx_strong,2);
datapool.median_onset_strong_highcond(find(datapool.median_onset_strong_highcond==0))= [];
datapool.median_onset_weak_highcond = datapool.median_onset(datapool.idx_weak,2);
datapool.median_onset_weak_highcond(find(datapool.median_onset_weak_highcond==0))= [];
datapool.median_onset_strong_highcond_highI = datapool.median_onset(datapool.idx_strong,4);
datapool.median_onset_strong_highcond_highI(find(datapool.median_onset_strong_highcond_highI==0))= [];
datapool.median_onset_weak_highcond_highI = datapool.median_onset(datapool.idx_weak,4);
datapool.median_onset_weak_highcond_highI(find(datapool.median_onset_weak_highcond_highI==0))= [];

datapool.sumarea_burst_lowcond = zeros(1, length(datapool.idx_burst_lowcond));
datapool.maxarea_burst_lowcond = zeros(1, length(datapool.idx_burst_lowcond));
datapool.sumarea_single_lowcond = zeros(1, length(datapool.idx_single_lowcond));
datapool.maxarea_single_lowcond = zeros(1, length(datapool.idx_single_lowcond));
for i = 1:length(datapool.idx_burst_lowcond)
    datapool.sumarea_burst_lowcond(i) = sum(data(datapool.idx_burst_lowcond(i)).area_single);
    datapool.maxarea_burst_lowcond(i) = max(data(datapool.idx_burst_lowcond(i)).area_single);
end
for i = 1:length(datapool.idx_single_lowcond)
    datapool.sumarea_single_lowcond(i) = sum(data(datapool.idx_single_lowcond(i)).area_single);
    datapool.maxarea_single_lowcond(i) = max(data(datapool.idx_single_lowcond(i)).area_single);
end

datapool.sumarea_burst_highcond = zeros(1, length(datapool.idx_burst_highcond));
datapool.maxarea_burst_highcond = zeros(1, length(datapool.idx_burst_highcond));
datapool.sumarea_single_highcond = zeros(1, length(datapool.idx_single_highcond));
datapool.maxarea_single_highcond = zeros(1, length(datapool.idx_single_highcond));
for i = 1:length(datapool.idx_burst_highcond)
    datapool.sumarea_burst_highcond(i) = sum(data(datapool.idx_burst_highcond(i)).area_single);
    datapool.maxarea_burst_highcond(i) = max(data(datapool.idx_burst_highcond(i)).area_single);
end
for i = 1:length(datapool.idx_single_highcond)
    datapool.sumarea_single_highcond(i) = sum(data(datapool.idx_single_highcond(i)).area_single);
    datapool.maxarea_single_highcond(i) = max(data(datapool.idx_single_highcond(i)).area_single);
end

datapool.sumarea_burst_highcond_highI = zeros(1, length(datapool.idx_burst_highcond_highI));
datapool.maxarea_burst_highcond_highI = zeros(1, length(datapool.idx_burst_highcond_highI));
datapool.sumarea_single_highcond_highI = zeros(1, length(datapool.idx_single_highcond_highI));
datapool.maxarea_single_highcond_highI = zeros(1, length(datapool.idx_single_highcond_highI));
for i = 1:length(datapool.idx_burst_highcond_highI)
    datapool.sumarea_burst_highcond_highI(i) = sum(data(datapool.idx_burst_highcond_highI(i)).area_single);
    datapool.maxarea_burst_highcond_highI(i) = max(data(datapool.idx_burst_highcond_highI(i)).area_single);
end
for i = 1:length(datapool.idx_single_highcond_highI)
    datapool.sumarea_single_highcond_highI(i) = sum(data(datapool.idx_single_highcond_highI(i)).area_single);
    datapool.maxarea_single_highcond_highI(i) = max(data(datapool.idx_single_highcond_highI(i)).area_single);
end

datapool.idx_highocnd_only = [];
for i = 1:length(data)
    if size(data(i).stim_spike_time,1)==1
        continue
    end
    clust = data(i).clust;
    clust(find(clust==999)) = [];
    if length(clust)<2
        continue
    end
    datapool.idx_highocnd_only = [datapool.idx_highocnd_only, i];
end
datapool.p_spike_highcond = datapool.p_spike(datapool.idx_highocnd_only,:);
datapool.p_burst_highcond = datapool.p_burst_true(datapool.idx_highocnd_only,:);
datapool.median_onset_highcond = datapool.median_onset(datapool.idx_highocnd_only,:);
%% pool together for bkg+clust
datapool.idx_4clust = [];
datapool.idx_3clust = [];
datapool.idx_2clust = [];
datapool.idx_1clust = [];

for i = 1:length(data)
    clust = data(i).clust;
    clust(find(clust==999)) = [];
    if length(data(i).isi)==1
        continue
    end
    if length(clust)==4
        datapool.idx_4clust = [datapool.idx_4clust, i];
    elseif length(clust)==3
        datapool.idx_3clust = [datapool.idx_3clust, i];
    elseif length(clust)==2
        datapool.idx_2clust = [datapool.idx_2clust, i];
    elseif length(clust)==1
        datapool.idx_1clust = [datapool.idx_1clust, i];
    end
end

datapool.pburst_4clust = datapool.p_burst(datapool.idx_4clust, 2);
datapool.pburst_3clust = datapool.p_burst(datapool.idx_3clust, 2);
datapool.pburst_2clust = datapool.p_burst(datapool.idx_2clust, 2);
datapool.pburst_1clust = datapool.p_burst(datapool.idx_1clust, 2);
datapool.APnum_4clust = datapool.AP_num(datapool.idx_4clust, 2);
datapool.APnum_3clust = datapool.AP_num(datapool.idx_3clust, 2);
datapool.APnum_2clust = datapool.AP_num(datapool.idx_2clust, 2);
datapool.APnum_1clust = datapool.AP_num(datapool.idx_1clust, 2);
datapool.isi_4clust = cell2mat(datapool.isi_highcond(datapool.idx_4clust));
datapool.min_isi_3clust = [];
for i = 1:length(datapool.idx_3clust)
    datapool.min_isi_3clust = [datapool.min_isi_3clust; min(datapool.isi_highcond{datapool.idx_3clust(i)})];
end
datapool.min_isi_2clust = [];
for i = 1:length(datapool.idx_2clust)
    datapool.min_isi_2clust = [datapool.min_isi_2clust; min(datapool.isi_highcond{datapool.idx_2clust(i)})];
end
datapool.min_isi_1clust = [];
for i = 1:length(datapool.idx_1clust)
    datapool.min_isi_1clust = [datapool.min_isi_1clust; min(datapool.isi_highcond{datapool.idx_1clust(i)})];
end
datapool.jitter_4clust = datapool.jitter(datapool.idx_4clust, 2);
datapool.jitter_3clust = datapool.jitter(datapool.idx_3clust, 2);
datapool.jitter_2clust = datapool.jitter(datapool.idx_2clust, 2);
datapool.jitter_1clust = datapool.jitter(datapool.idx_1clust, 2);

figure
g1 = repmat({'a'},length(datapool.idx_1clust),1);
g2 = repmat({'b'},length(datapool.idx_2clust),1);
g3 = repmat({'c'},length(datapool.idx_3clust),1);
y = [datapool.pburst_1clust; datapool.pburst_2clust; datapool.pburst_3clust];
boxplot(y, [g1;g2;g3], 'Color', ['r','b' ,'k'],'PlotStyle','compact')
%% jitter of the second spikes evoked by patterns containing a Na d-spike
idx = intersect([datapool.idx_4clust, datapool.idx_3clust, datapool.idx_2clust], datapool.idx_strong);
jitter_second_spike = [];
data_temp = data(idx);
first_jitter = [];
second_jitter = [];
patterns = [];
for i  = 1:length(data_temp)
    for n = 2:5
        if length(data_temp(i).p_burst)>=n
            if data_temp(i).p_burst(n)>=0.2
                patterns = [patterns,i];
                second_spike = [];
                first_spike = [];
                for j = 1:10
                    if length(data_temp(i).stim_spike_time{n,j})>=1
                        first_spike = [first_spike, data_temp(i).stim_spike_time{n,j}(1)];
                    end
                    if length(data_temp(i).stim_spike_time{n,j})>=2
                        second_spike = [second_spike, data_temp(i).stim_spike_time{n,j}(2)];
                    end
                end
                first_jitter = [first_jitter,std(first_spike)];
                second_jitter = [second_jitter,std(second_spike)];
            end
        end
    end
    
end
%%
datapool.pburst_4clust_lowcond = [];
datapool.APnum_4clust_lowcond = [];
datapool.isi_4clust_lowcond = [];
datapool.onset_4clust_lowcond = [];

datapool.pburst_3clust_lowcond = [];
datapool.APnum_3clust_lowcond = [];
datapool.isi_3clust_lowcond = [];
datapool.onset_3clust_lowcond = [];

datapool.pburst_2clust_lowcond = [];
datapool.APnum_2clust_lowcond = [];
datapool.isi_2clust_lowcond = [];
datapool.onset_2clust_lowcond = [];

datapool.pburst_1clust_lowcond = [];
datapool.APnum_1clust_lowcond = [];
datapool.isi_1clust_lowcond = [];
datapool.onset_1clust_lowcond = [];
idx = [];
for i = 1:length(datapool.samecell_idx)
    
    idx_1clust = [];
    idx_2clust = [];
    idx_3clust = [];
    for j = 1:length(datapool.samecell_idx{i})

        clust = data(datapool.samecell_idx{i}(j)).clust;
        clust(find(clust==999)) = [];
        if (length(clust)==1) 
            idx_1clust = [idx_1clust,datapool.samecell_idx{i}(j)];
        elseif (length(clust)==2) && (data(datapool.samecell_idx{i}(j)).p_spike(1)>0)
            idx_2clust = [idx_2clust,datapool.samecell_idx{i}(j)];
        elseif (length(clust)==3) && (data(datapool.samecell_idx{i}(j)).p_spike(1)>0)
            idx_3clust = [idx_3clust,datapool.samecell_idx{i}(j)];
        end
    end

    if (isempty(idx_3clust))||(isempty(idx_2clust))||(isempty(idx_1clust))
        continue
    end
    for j = idx_1clust
        clust_1 = data(j).clust;
        clust_1(find(clust_1==999)) = [];
        idx_2 = [];
        idx_3 = [];
        for n = idx_2clust
            clust_2 = data(n).clust;
            clust_2(find(clust_2==999)) = [];
            if isempty(find(clust_2==clust_1))
                continue
            else
                idx_2 = [idx_2,n];
            end
        end
        for n = idx_3clust
            clust_3 = data(n).clust;
            clust_3(find(clust_3==999)) = [];
            if isempty(find(clust_3==clust_1))
                continue
            else
                idx_3 = [idx_3,n];
            end
        end
        if (isempty(idx_2))||(isempty(idx_3))
            continue
        end
        if length(idx_2)>1
            order = [];
            for n = idx_2
                clust_2 = data(n).clust;
                clust_2(find(clust_2==999)) = [];
                order = [order, find(clust_2==clust_1)];
            end
            [~,a] = min(idx_2);
            idx_2 = idx_2(a);
        end
        if length(idx_3)>1
            order = [];
            for n = idx_3
                clust_3 = data(n).clust;
                clust_3(find(clust_3==999)) = [];
                order = [order, find(clust_3==clust_1)];
            end
            [~,a] = min(idx_3);
            idx_3 = idx_3(a);
        end
        datapool.pburst_1clust_lowcond = [datapool.pburst_1clust_lowcond, data(j).p_burst(1)];
        datapool.APnum_1clust_lowcond = [datapool.APnum_1clust_lowcond, data(j).AP_num(1)];
        datapool.isi_1clust_lowcond = [datapool.isi_1clust_lowcond, min(data(j).isi{1})];
        datapool.onset_1clust_lowcond = [datapool.onset_1clust_lowcond, data(j).onset(1)];
        datapool.pburst_2clust_lowcond = [datapool.pburst_2clust_lowcond, data(idx_2).p_burst(1)];
        datapool.APnum_2clust_lowcond = [datapool.APnum_2clust_lowcond, data(idx_2).AP_num(1)];
        datapool.isi_2clust_lowcond = [datapool.isi_2clust_lowcond, min(data(idx_2).isi{1})];
        datapool.onset_2clust_lowcond = [datapool.onset_2clust_lowcond, data(idx_2).onset(1)];
        datapool.pburst_3clust_lowcond = [datapool.pburst_3clust_lowcond, data(idx_3).p_burst(1)];
        datapool.APnum_3clust_lowcond = [datapool.APnum_3clust_lowcond, data(idx_3).AP_num(1)];
        datapool.isi_3clust_lowcond = [datapool.isi_3clust_lowcond, min(data(idx_3).isi{1})];
        datapool.onset_3clust_lowcond = [datapool.onset_3clust_lowcond, data(idx_3).onset(1)];
        idx = [idx, i];
    end
end
color1 = [241,90,36]/256;
color2 = [0,164,79]/255;
color3 = [0,113,188]/256;
figure
y = [datapool.pburst_1clust_lowcond', datapool.pburst_2clust_lowcond',datapool.pburst_3clust_lowcond'];
boxplot(y, 'Color', [color3;color2;color1],'PlotStyle','compact')
%% decode - sort dataset

% sort dataset
[x,idx]=sortrows([{data(datarange).date}; {data(datarange).cell}]', [1,2]);
data1 = data(datarange);
data1 = data1(idx);
if datarange(1)==1
    data = data1;
else
    data = [data(1:datarange(1)-1), data1];
end
%% decode - build stim responser kernel
datarange = 1:length(data);
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
%% detect same cells, build decode dataset
if datarange(1) ==1
    datapool.samecell_idx = {[1]};
    count = 1;
    for i = 2:length(data)
        if (str2num(data(i).date)~=str2num(data(i-1).date))|(data(i).cell~=data(i-1).cell)
            count = count + 1;
            datapool.samecell_idx{count} = [i];
        else
            datapool.samecell_idx{count} = [datapool.samecell_idx{count}, i];
        end
    end
    decoderange = 1:length(datapool.samecell_idx);
else
    decoderange = length(datapool.samecell_idx) +1;
    datapool.samecell_idx = [datapool.samecell_idx, {[datarange(1)]}];
    count = length(datapool.samecell_idx);
    for i = datarange(2:end)
        if (str2num(data(i).date)~=str2num(data(i-1).date))|(data(i).cell~=data(i-1).cell)
            count = count + 1;
            datapool.samecell_idx{count} = [i];
        else
            datapool.samecell_idx{count} = [datapool.samecell_idx{count}, i];
        end
    end
    decoderange = decoderange:length(datapool.samecell_idx);
end


%% build kernel for spontansoue spikes    
for i = decoderange
    data_decode(i).date = data(datapool.samecell_idx{i}(1)).date;
    data_decode(i).cell = data(datapool.samecell_idx{i}(1)).cell;
    data_decode(i).idx = datapool.samecell_idx{i};
    data_decode(i).clust = {data(datapool.samecell_idx{i}).clust};
    data_decode(i).strong_clust = [];
    for j = 1:length(datapool.samecell_idx{i})
        idx_strong = find((data(datapool.samecell_idx{i}(j)).maxdv_single>2)&((data(datapool.samecell_idx{i}(j)).maxdv_single<50)));
        data_decode(i).strong_num(j) = length(idx_strong);
        if data_decode(i).strong_num(j) >= 1
            data_decode(i).strong_clust = [data_decode(i).strong_clust,data(datapool.samecell_idx{i}(j)).clust(idx_strong)];
        end        
    end
    data_decode(i).strong_clust = unique(data_decode(i).strong_clust);
    for j = 1:5
        data_decode(i).stim_spike_time_matrix{j} = [];
        data_decode(i).stim_spike_time_kernel{j} = [];
        for n = 1:length(datapool.samecell_idx{i})
            data_decode(i).stim_spike_time_matrix{j} = [data_decode(i).stim_spike_time_matrix{j}; data(datapool.samecell_idx{i}(n)).stim_spike_time_matrix{j}];
            data_decode(i).stim_spike_time_kernel{j} = [data_decode(i).stim_spike_time_kernel{j}; data(datapool.samecell_idx{i}(n)).stim_spike_time_kernel{j}];
        end
    end
    for j = 1:5
        data_decode(i).spon_spike_time{j} = [];
        trial_break = 0;
        for n = 1:length(datapool.samecell_idx{i})
            if size(data(datapool.samecell_idx{i}(n)).stim_spike_time, 1)<j
                continue
            else
                spike_time = data(datapool.samecell_idx{i}(n)).spike_time{j};
                orig_spike_time = data(datapool.samecell_idx{i}(n)).spike_time{j};
                flag = 0;
                stimtrig_idx = [];
                for m = 1:length(data(datapool.samecell_idx{i}(n)).stim_start_time{j})
                    stimtrig_idx = [stimtrig_idx, find((orig_spike_time-data(datapool.samecell_idx{i}(n)).stim_start_time{j}(m)<=win_length)&(orig_spike_time-data(datapool.samecell_idx{i}(n)).stim_start_time{j}(m)>0))];
                    if ~isempty(find((orig_spike_time-data(datapool.samecell_idx{i}(n)).stim_start_time{j}(m)<=win_length)&(orig_spike_time-data(datapool.samecell_idx{i}(n)).stim_start_time{j}(m)>0)))
                        spike_time(stimtrig_idx(end)+1:end) = spike_time(stimtrig_idx(end)+1:end)-win_length;
                        flag = flag+1;
                    end
                end
                spike_time(stimtrig_idx) = [];
                data_decode(i).spon_spike_time{j} = [data_decode(i).spon_spike_time{j}, spike_time + trial_break];
                trial_break = trial_break + 10-win_length*flag;
            end
            
        end
        if isempty(data_decode(i).spon_spike_time{j})
            data_decode(i).spon_spike_time_matrix{j} = zeros(1,1e4);
            data_decode(i).spon_spike_time_kernel{j} = zeros(1,1e4);
            continue
        end
        edges = 0:bin_size:ceil(data_decode(i).spon_spike_time{j}(end)/win_length)*win_length;
        data_decode(i).spon_spike_time_matrix{j} = histcounts(data_decode(i).spon_spike_time{j}, edges);
        if length(data_decode(i).spon_spike_time_matrix{j}) < 1e4
            data_decode(i).spon_spike_time_matrix{j} = [data_decode(i).spon_spike_time_matrix{j}, zeros(1, 1e4-length(data_decode(i).spon_spike_time_matrix{j}))];
        end
        a = conv(data_decode(i).spon_spike_time_matrix{j},w);
        data_decode(i).spon_spike_time_kernel{j} = a(1:length(data_decode(i).spon_spike_time_matrix{j}));
        
    end    
end
for i = decoderange
    data_decode(i).dataset = cell(1,5);
    for j = 1:length(data_decode(i).idx)
        N = size(data(data_decode(i).idx(j)).stim_spike_time, 1);
        for n = 1:N
            data_decode(i).dataset{n} = [data_decode(i).dataset{n}, j];
        end
    end
end
for i = decoderange
    data_decode(i).spon_FR = [];
    for j = 1:length(data_decode(i).spon_spike_time_matrix)
        data_decode(i).spon_FR(j) = length(find(data_decode(i).spon_spike_time_matrix{j}~=0))/(length(data_decode(i).spon_spike_time_matrix{j})*bin_size);
    end
    datapool.spon_FR(i,:) = data_decode(i).spon_FR;
end
%% cross correlation
datapool.n_trials = 10;
for i = decoderange
    data_decode(i).correlation_matrix = cell(1,5);
    data_decode(i).correlation_ave = cell(1,5);
    for j = 1:5
        a = [];
        if isempty(data_decode(i).dataset{j})
            continue
        end
        for n = 1:length(data_decode(i).dataset{j})
            a = [a, (data_decode(i).dataset{j}(n)-1)*datapool.n_trials+1:data_decode(i).dataset{j}(n)*datapool.n_trials];
        end
        A = data_decode(i).stim_spike_time_kernel{j}(a,:);
        A = A';
%         A = A-mean(A);
        A1 = corr(A);
        A1(isnan(A1))=0;
        A_sort = cov(A);
        A2 = zeros(size(A));
        N = round(size(A1, 1)/datapool.n_trials);
        data_decode(i).correlation_ave{j} = zeros(N,N);
        for n = 1:N
            A_temp = A(:,(n-1)*datapool.n_trials+1:n*datapool.n_trials);
            corr_temp = cov(A_temp);
            I = 1:datapool.n_trials;
            for m = 1:datapool.n_trials
                [~,a] = sort(corr_temp(m:end, m), 'descend');
                I(m:end) = a+m-1;
                A_temp(:,m:end) = A_temp(:,a+m-1);
                corr_temp = cov(A_temp);
            end
            A2(:,(n-1)*datapool.n_trials+1:n*datapool.n_trials) = A_temp;
            data_decode(i).correlation_matrix{j} = corr(A2);
            for m = 1:N
                if m == n
                    data_decode(i).correlation_ave{j}(n,m) = (sum(sum(A1((n-1)*datapool.n_trials+1:n*datapool.n_trials, (n-1)*datapool.n_trials+1:n*datapool.n_trials)))-sum(diag(A1((n-1)*datapool.n_trials+1:n*datapool.n_trials, (n-1)*datapool.n_trials+1:n*datapool.n_trials))))/(datapool.n_trials^2-datapool.n_trials);
                    data_decode(i).correlation_max{j}(n,m) = max(max(triu(A1((n-1)*datapool.n_trials+1:n*datapool.n_trials, (n-1)*datapool.n_trials+1:n*datapool.n_trials), 1)));
                else
                    data_decode(i).correlation_ave{j}(n,m) = sum(sum(A1((n-1)*datapool.n_trials+1:n*datapool.n_trials, (m-1)*datapool.n_trials+1:m*datapool.n_trials)))/(datapool.n_trials^2);
                    data_decode(i).correlation_max{j}(n,m) = max(max(A1((n-1)*datapool.n_trials+1:n*datapool.n_trials, (m-1)*datapool.n_trials+1:m*datapool.n_trials)));
                end
            end
        end        

    end 
%         figure
%         imagesc(correlation_matrix{j,i})
%         colormap('jet')
%         caxis([-0.1,0.1])
end

%% pool of cross correlation
datapool.corr_pool = cell(5,1);
datapool.corr_auto = cell(5,1);
datapool.corr_cross = cell(5,1);
datapool.corr_pool_byclust = cell(5,1);
for j = 1:5
    datapool.corr_pool{j} = zeros(length(datapool.samecell_idx), 2);
    datapool.corr_pool_byclust{j} = [];
    datapool.corr_auto{j} = [];
    datapool.corr_cross{j} = [];
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{j})
            continue
        end
        datapool.corr_pool{j}(i,1) = max(diag(data_decode(i).correlation_ave{j}));
        temp_matrix = data_decode(i).correlation_ave{j};
        for n = 1:size(temp_matrix, 1)
            temp_matrix(n,n) = 0;
        end
        datapool.corr_pool{j}(i,2) = max(max(temp_matrix));
        datapool.corr_pool_byclust{j} = [datapool.corr_pool_byclust{j};[diag(data_decode(i).correlation_ave{j}),sum(temp_matrix,2)/(datapool.n_trials-1)]];
        temp_matrix = [];
        for n = 1:length(data_decode(i).dataset{j})
            for m = n: length(data_decode(i).dataset{j})
                A_temp =  data_decode(i).correlation_matrix{j}((n-1)*datapool.n_trials+1:n*datapool.n_trials, (m-1)*datapool.n_trials+1:m*datapool.n_trials);
                if m==n
                   a = reshape(triu(A_temp, 1), 1, datapool.n_trials^2); 
                   a(find(a==0)) = [];
                   datapool.corr_auto{j} = [datapool.corr_auto{j}, a];
                else
                   a = reshape(triu(A_temp), 1, datapool.n_trials^2);
                   a(find(a==0)) = [];
                   datapool.corr_cross{j} = [datapool.corr_cross{j}, a];
                end
            end
        end
    end
end


%% discrimination between clusters
for i = decoderange
    data_decode(i).linear_regression_R = cell(1,5);
    data_decode(i).linear_regression_error = cell(1, 5);
    data_decode(i).linear_regression_R_shuffle_clust = cell(1, 5);
    data_decode(i).linear_regression_error_shuffle_clust = cell(1, 5);
    for j = 1:5     
        a = [];
        if isempty(data_decode(i).dataset{j})
            continue
        end
        for n = 1:length(data_decode(i).dataset{j})
            a = [a, (data_decode(i).dataset{j}(n)-1)*datapool.n_trials+1:data_decode(i).dataset{j}(n)*datapool.n_trials];
        end
        A = data_decode(i).stim_spike_time_kernel{j}(a,:);
        [data_decode(i).linear_regression_R{j}, data_decode(i).linear_regression_error{j}] = multicluster_linear_regression(A, datapool.n_trials);
        [data_decode(i).linear_regression_R_shuffle_clust{j}, data_decode(i).linear_regression_error_shuffle_clust{j}] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
    for n = 1:9
        [R_shuffle_clust, error_shuffle_clust] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
        data_decode(i).linear_regression_R_shuffle_clust{j} = data_decode(i).linear_regression_R_shuffle_clust{j} + R_shuffle_clust;
        data_decode(i).linear_regression_error_shuffle_clust{j} = data_decode(i).linear_regression_error_shuffle_clust{j} + error_shuffle_clust;
    end
    data_decode(i).linear_regression_R_shuffle_clust{j} = data_decode(i).linear_regression_R_shuffle_clust{j}/10;
    data_decode(i).linear_regression_error_shuffle_clust{j} = data_decode(i).linear_regression_error_shuffle_clust{j}/10;
    end
    disp(sprintf('cell number %d done', i));
end

%% discrimination between cluster and spontaneous spike
for i = decoderange
    data_decode(i).linear_regression_R_stim_spon = cell(size(data_decode(i).clust,2),5);
    data_decode(i).linear_regression_R_stim_spon_shuffle = cell(size(data_decode(i).clust,2),5);
    for j = 1:5     
        if isempty(data_decode(i).dataset{j})
            continue
        end
        count = 0;
        A_spon = [];
        while (size(A_spon,1)<20)&&count<10000
            m = randsample(length(data_decode(i).spon_spike_time_kernel{j})-ceil(win_length/bin_size),1);
            if max(data_decode(i).spon_spike_time_kernel{j}(m:m+ceil(win_length/bin_size)-1))>0
                A_spon = [A_spon;data_decode(i).spon_spike_time_kernel{j}(m:m+ceil(win_length/bin_size)-1)];
            end
            count = count + 1;
        end 
        if size(A_spon, 1)<20
            A_spon = [A_spon;zeros(20-size(A_spon,1),ceil(win_length/bin_size))];
        end
        A_spon = A_spon(randperm(size(A_spon, 1)),:);
        data_decode(i).spon_trials{j} = A_spon;
        for n = 1:length(data_decode(i).dataset{j})
            A = data_decode(i).stim_spike_time_kernel{j}((data_decode(i).dataset{j}(n)-1)*datapool.n_trials+1:data_decode(i).dataset{j}(n)*datapool.n_trials,:);
            A = [A;A_spon];
            [data_decode(i).linear_regression_R_stim_spon{n,j},~] = multicluster_linear_regression(A, datapool.n_trials);
            [data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}, ~] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
            for m = 1:9
                [R_shuffle, ~] = multicluster_linear_regression_shuffle(A, datapool.n_trials);
                data_decode(i).linear_regression_R_stim_spon_shuffle{n,j} = data_decode(i).linear_regression_R_stim_spon_shuffle{n,j} + R_shuffle;
            end
            data_decode(i).linear_regression_R_stim_spon_shuffle{n,j} = data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}/10;
        end
    end
    disp(sprintf('cell number %d done', i));
end
%% pool of discriminability
datapool.R_pool = cell(5,1);
for n = 1:5
    datapool.R_pool{n} = zeros(length(datapool.samecell_idx), 2);
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{n})
            continue
        end
        datapool.R_pool{n}(i,1) = max(max(data_decode(i).linear_regression_R{n}));
        datapool.R_pool{n}(i,2) = max(max(data_decode(i).linear_regression_R_shuffle_clust{n}));
    end
end
datapool.error_pool = cell(5,1);
for n = 1:5
        datapool.error_pool{n} = zeros(length(datapool.samecell_idx), 2);
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{n})
            continue
        end
        datapool.error_pool{n}(i,1) = max(max(data_decode(i).linear_regression_error{n}));
        datapool.error_pool{n}(i,2) = max(max(data_decode(i).linear_regression_error_shuffle_clust{n}));
    end
end
[h,p] = signrank(datapool.R_pool{1}(:,1), datapool.R_pool{1}(:,2))
[h,p] = signrank([datapool.R_pool{2}(:,1); datapool.R_pool{3}(:,1)], [datapool.R_pool{2}(:,2); datapool.R_pool{3}(:,2)])
[h,p] = signrank([datapool.R_pool{4}(:,1); datapool.R_pool{5}(:,1)], [datapool.R_pool{4}(:,2); datapool.R_pool{5}(:,2)])

%% R between clusters as a function of firing rate
bin_size = 2;
freq_bin = [0:bin_size:10];
freq_c = bin_size/2+freq_bin(1:end-1);
datapool.R_FR = cell(1,length(freq_c));
datapool.error_FR = cell(1,length(freq_c));
for i = 1:length(datapool.samecell_idx)
    for j = 2:5
        for n = 1:length(data_decode(i).clust)-1
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
                if size(data_decode(i).linear_regression_R{j}, 1)<max(n,m)
                    continue
                end
                bkg_FR_1 = length(data(data_decode(i).idx(n)).spike_time{j})/10;
                bkg_FR_2 = length(data(data_decode(i).idx(m)).spike_time{j})/10;
                bkg_FR_ave = mean([bkg_FR_1, bkg_FR_2]);
                [temp, ~] = histcounts(bkg_FR_ave, freq_bin);
                if isempty(find(temp>0))
                    continue
                else
                    datapool.R_FR{find(temp>0)} = [datapool.R_FR{find(temp>0)};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                    datapool.error_FR{find(temp>0)} = [datapool.error_FR{find(temp>0)};[data_decode(i).linear_regression_error{j}(n,m), data_decode(i).linear_regression_error_shuffle_clust{j}(n,m)]];
                end
                    
            end
        end

    end
end
g = [];
R_temp = [];
R_temp_s = [];
for i = 1:length(freq_c)
    g = [g; repmat(freq_c(i),size(datapool.R_FR{i},1),1)];
    R_temp = [R_temp;datapool.R_FR{i}(:,1)];
    R_temp_s = [R_temp_s;datapool.R_FR{i}(:,2)];
end
figure
boxplot(R_temp, g, 'PlotStyle','compact','Color',[240, 90, 36]/255)
ylim([0,1])
figure
boxplot(R_temp_s, g, 'PlotStyle','compact','Color',[128, 128, 128]/255)
ylim([0,1])


R_mean = zeros(2, length(freq_c));
R_sem = zeros(2, length(freq_c));
error_mean = zeros(2, length(freq_c));
error_sem = zeros(2, length(freq_c));
for i = 1:length(freq_c)
    R_mean(:,i) = mean(datapool.R_FR{i})';
    R_sem(:,i) = std(datapool.R_FR{i},[],1)'/size(datapool.R_FR{i},1);
    error_mean(:,i) = mean(datapool.error_FR{i})';
    error_sem(:,i) = std(datapool.error_FR{i},[],1)'/size(datapool.R_FR{i},1);
end

figure
errorbar(freq_c, R_mean(1,:), R_sem(1,:), 'o', 'MarkerSize',10,  'Color', [240, 90, 36]/255,'MarkerEdgeColor', [240, 90, 36]/255,'MarkerFaceColor', [1,1,1] )
hold on
errorbar(freq_c, R_mean(2,:), R_sem(2,:), 'o', 'MarkerSize',10,  'Color', [0.5,0.5,0.5],'MarkerEdgeColor', [0.5,0.5,0.5],'MarkerFaceColor', [1,1,1] )
hold on

figure
plot(freq_c, R_mean(1,:), 'Color', [240, 90, 36]/255,'Linewidth', 1)
hold on
plot(freq_c, R_mean(1,:)+R_sem(1,:), 'Color', [240, 90, 36]/255,'Linewidth', 0.5)
plot(freq_c, R_mean(1,:)-R_sem(1,:), 'Color', [240, 90, 36]/255,'Linewidth', 0.5)
plot(freq_c, R_mean(2,:), 'Color', [128, 128, 128]/255,'Linewidth', 1)
hold on
plot(freq_c, R_mean(2,:)+R_sem(2,:), 'Color', [128, 128, 128]/255,'Linewidth', 0.5)
plot(freq_c, R_mean(2,:)-R_sem(2,:), 'Color', [128, 128, 128]/255,'Linewidth', 0.5)
xlim([0,10])
ylim([0,1])
% 
% figure
% errorbar(freq_c, error_mean(1,:), error_sem(1,:), 'o', 'MarkerSize',10,  'Color', 'r','MarkerEdgeColor', 'r','MarkerFaceColor', [1,1,1] )
% hold on
% errorbar(freq_c, error_mean(2,:), error_sem(2,:), 'o', 'MarkerSize',10,  'Color', [0.5,0.5,0.5],'MarkerEdgeColor', [0.5,0.5,0.5],'MarkerFaceColor', [1,1,1] )
% hold on
%% discrimination between strong and weak pattern
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

y1 = cell2mat(datapool.R_strong(2:5)');
y2 = cell2mat(datapool.R_strongweak(2:5)');


g1 = repmat({'a'},size(y1,1),1);
g11 = repmat({'a1'},size(y1,1),1);
g2 = repmat({'b'},size(y2,1),1);
g21 = repmat({'b1'},size(y2,1),1);

figure
boxplot([y1(:,1);y1(:,2);y2(:,1);y2(:,2)], [g1;g11;g2;g21], 'Color', [datapool.colorpool(:,:,1);[0.5,0.5,0.5];datapool.colorpool(:,:,2); [0.5,0.5,0.5];],'PlotStyle','compact')
%% discrimination between strong and weak pattern (4 clust only)
datapool.R_strong_4clust = cell(1,5);
datapool.R_weak_4clust = cell(1,5);
datapool.R_strongweak_4clust = cell(1,5);
for i = 1:length(datapool.samecell_idx)
    for j = 1:5
        for n = 1:length(data_decode(i).clust)-1
            if size(data_decode(i).linear_regression_R{j}, 1)<n
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)<data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(1,2)
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)==0
                continue
            end
            clust = data_decode(i).clust{n};
            clust(find(clust==999)) = [];
            if length(clust)<4
                continue;
            end
            for m = n+1:length(data_decode(i).clust)
                if size(data_decode(i).linear_regression_R{j}, 1)<max(n,m)
                    continue
                end
                clust = data_decode(i).clust{m};
                clust(find(clust==999)) = [];
                if length(clust)<4
                    continue;
                end
                if ((data_decode(i).strong_num(n)>=1)&&(data_decode(i).strong_num(m)<1))||((data_decode(i).strong_num(n)<1)&&(data_decode(i).strong_num(m)>=1))
                    datapool.R_strongweak_4clust{j} = [datapool.R_strongweak_4clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((data_decode(i).strong_num(n)>=1)&&(data_decode(i).strong_num(m)>=1))
                     datapool.R_strong_4clust{j} = [datapool.R_strong_4clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((data_decode(i).strong_num(n)<1)&&(data_decode(i).strong_num(m)<1))
                    datapool.R_weak_4clust{j} = [datapool.R_weak_4clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                end
            end
        end

    end
end
g1 = repmat({'a'},size(datapool.R_strong_4clust{2},1)+size(datapool.R_strong_4clust{3},1)+size(datapool.R_strong_4clust{4},1)+size(datapool.R_strong_4clust{5},1),1);
g11 = repmat({'a1'},size(datapool.R_strong_4clust{2},1)+size(datapool.R_strong_4clust{3},1)+size(datapool.R_strong_4clust{4},1)+size(datapool.R_strong_4clust{5},1),1);
g2 = repmat({'b'},size(datapool.R_weak_4clust{2},1)+size(datapool.R_weak_4clust{3},1)+size(datapool.R_weak_4clust{4},1)+size(datapool.R_weak_4clust{5},1),1);
g21 = repmat({'b1'},size(datapool.R_weak_4clust{2},1)+size(datapool.R_weak_4clust{3},1)+size(datapool.R_weak_4clust{4},1)+size(datapool.R_weak_4clust{5},1),1);
g3 = repmat({'c'},size(datapool.R_strongweak_4clust{2},1)+size(datapool.R_strongweak_4clust{3},1)+size(datapool.R_strongweak_4clust{4},1)+size(datapool.R_strongweak_4clust{5},1),1);
g31 = repmat({'c1'},size(datapool.R_strongweak_4clust{2},1)+size(datapool.R_strongweak_4clust{3},1)+size(datapool.R_strongweak_4clust{4},1)+size(datapool.R_strongweak_4clust{5},1),1);
figure
boxplot([datapool.R_strong_4clust{2}(:,1);datapool.R_strong_4clust{3}(:,1);datapool.R_strong_4clust{4}(:,1);datapool.R_strong_4clust{5}(:,1);
    datapool.R_strong_4clust{2}(:,2);datapool.R_strong_4clust{3}(:,2);datapool.R_strong_4clust{4}(:,2);datapool.R_strong_4clust{5}(:,2);
    datapool.R_weak_4clust{2}(:,1);datapool.R_weak_4clust{3}(:,1);datapool.R_weak_4clust{4}(:,1);datapool.R_weak_4clust{5}(:,1);
    datapool.R_weak_4clust{2}(:,2);datapool.R_weak_4clust{3}(:,2);datapool.R_weak_4clust{4}(:,2);datapool.R_weak_4clust{5}(:,2);
    datapool.R_strongweak_4clust{2}(:,1);datapool.R_strongweak_4clust{3}(:,1);datapool.R_strongweak_4clust{4}(:,1);datapool.R_strongweak_4clust{5}(:,1);
    datapool.R_strongweak_4clust{2}(:,2);datapool.R_strongweak_4clust{3}(:,2);datapool.R_strongweak_4clust{4}(:,2);datapool.R_strongweak_4clust{5}(:,2)], [g1;g11;g2;g21;g3;g31],  'Color', [datapool.colorpool(:,:,1);[0,0,0];datapool.colorpool(:,:,2);[0,0,0];datapool.colorpool(:,:,3);[0,0,0]],'PlotStyle','compact')


% g1 = repmat({'a'},size(datapool.R_strong_4clust{2},1)+size(datapool.R_strong_4clust{3},1),1);
% g11 = repmat({'a1'},size(datapool.R_strong_4clust{2},1)+size(datapool.R_strong_4clust{3},1),1);
% g2 = repmat({'b'},size(datapool.R_weak_4clust{2},1)+size(datapool.R_weak_4clust{3},1),1);
% g21 = repmat({'b1'},size(datapool.R_weak_4clust{2},1)+size(datapool.R_weak_4clust{3},1),1);
% g3 = repmat({'c'},size(datapool.R_strongweak_4clust{2},1)+size(datapool.R_strongweak_4clust{3},1),1);
% g31 = repmat({'c1'},size(datapool.R_strongweak_4clust{2},1)+size(datapool.R_strongweak_4clust{3},1),1);
% figure
% boxplot([datapool.R_strong_4clust{2}(:,1);datapool.R_strong_4clust{3}(:,1);
%     datapool.R_strong_4clust{2}(:,2);datapool.R_strong_4clust{3}(:,2);
%     datapool.R_weak_4clust{2}(:,1);datapool.R_weak_4clust{3}(:,1);
%     datapool.R_weak_4clust{2}(:,2);datapool.R_weak_4clust{3}(:,2);
%     datapool.R_strongweak_4clust{2}(:,1);datapool.R_strongweak_4clust{3}(:,1);
%     datapool.R_strongweak_4clust{2}(:,2);datapool.R_strongweak_4clust{3}(:,2)], [g1;g11;g2;g21;g3;g31],  'Color', [datapool.colorpool(:,:,1);[0,0,0];datapool.colorpool(:,:,2);[0,0,0];datapool.colorpool(:,:,3);[0,0,0]],'PlotStyle','compact')

% g1 = repmat({'a'},size(datapool.R_strong_4clust{2},1),1);
% g11 = repmat({'a1'},size(datapool.R_strong_4clust{2},1),1);
% g2 = repmat({'b'},size(datapool.R_weak_4clust{2},1),1);
% g21 = repmat({'b1'},size(datapool.R_weak_4clust{2},1),1);
% g3 = repmat({'c'},size(datapool.R_strongweak_4clust{2},1),1);
% g31 = repmat({'c1'},size(datapool.R_strongweak_4clust{2},1),1);
% figure
% boxplot([datapool.R_strong_4clust{2}(:,1);
%     datapool.R_strong_4clust{2}(:,2);
%     datapool.R_weak_4clust{2}(:,1);
%     datapool.R_weak_4clust{2}(:,2);
%     datapool.R_strongweak_4clust{2}(:,1);
%     datapool.R_strongweak_4clust{2}(:,2)], [g1;g11;g2;g21;g3;g31],  'Color', [datapool.colorpool(:,:,1);[0,0,0];datapool.colorpool(:,:,2);[0,0,0];datapool.colorpool(:,:,3);[0,0,0]],'PlotStyle','compact')

%% discrimination between patterns with different number of clusters
datapool.R_4clust = cell(1,5);
datapool.R_4_3clust = cell(1,5);
datapool.R_2_3clust = cell(1,5);
datapool.R_1_2clust = cell(1,5);
datapool.R_1_3clust = cell(1,5);
for i = 1:length(datapool.samecell_idx)
    for j = 1:5
        for n = 1:length(data_decode(i).clust)-1
            if size(data_decode(i).linear_regression_R{j}, 1)<n
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)<data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(1,2)
                continue
            end
            if data_decode(i).linear_regression_R_stim_spon{n,j}(1,2)==0
                continue
            end
            nclust = data_decode(i).clust{n};
            nclust(find(nclust==999)) = [];
            nclust = length(nclust);
            for m = n+1:length(data_decode(i).clust)
                if size(data_decode(i).linear_regression_R{j}, 1)<max(n,m)
                    continue
                end
                mclust = data_decode(i).clust{m};
                mclust(find(mclust==999)) = [];
                mclust = length(mclust);
                if (mclust==4) && (nclust==4)
                    datapool.R_4clust{j} = [datapool.R_4clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((mclust==3)&&(nclust==4))||((mclust==4)&&(nclust==3))
                    datapool.R_4_3clust{j} = [datapool.R_4_3clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((mclust==3)&&(nclust==2))||((mclust==2)&&(nclust==3))
                    datapool.R_2_3clust{j} = [datapool.R_2_3clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((mclust==2)&&(nclust==1))||((mclust==1)&&(nclust==2))
                    datapool.R_1_2clust{j} = [datapool.R_1_2clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                elseif ((mclust==3)&&(nclust==1))||((mclust==1)&&(nclust==3))
                    datapool.R_1_3clust{j} = [datapool.R_1_3clust{j};[data_decode(i).linear_regression_R{j}(n,m), data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
                end
            end
        end

    end
end


g1 = repmat({'a'},size(datapool.R_1_2clust{2},1)+size(datapool.R_1_2clust{3},1),1);
g11 = repmat({'a1'},size(datapool.R_1_2clust{2},1)+size(datapool.R_1_2clust{3},1),1);
g2 = repmat({'b'},size(datapool.R_1_3clust{2},1)+size(datapool.R_1_3clust{3},1),1);
g21 = repmat({'b1'},size(datapool.R_1_3clust{2},1)+size(datapool.R_1_3clust{3},1),1);
g3 = repmat({'c'},size(datapool.R_4_3clust{2},1)+size(datapool.R_4_3clust{3},1),1);
g31 = repmat({'c1'},size(datapool.R_4_3clust{2},1)+size(datapool.R_4_3clust{3},1),1);

y1 = [datapool.R_1_2clust{2}(:,1);datapool.R_1_2clust{3}(:,1)];
y11 = [datapool.R_1_2clust{2}(:,2);datapool.R_1_2clust{3}(:,2)];
y2 = [datapool.R_1_3clust{2}(:,1);datapool.R_1_3clust{3}(:,2)];
y21 = [datapool.R_1_3clust{2}(:,2);datapool.R_1_3clust{3}(:,2)];
y3 = [datapool.R_4_3clust{2}(:,1);datapool.R_4_3clust{3}(:,1)];
y31 = [datapool.R_4_3clust{2}(:,2);datapool.R_4_3clust{3}(:,2)];

figure
boxplot([y1;y11;y2;y21;y3;y31], [g1;g11;g2;g21;g3;g31],  'Color', [datapool.colorpool(:,:,1);[0,0,0];datapool.colorpool(:,:,2);[0,0,0];datapool.colorpool(:,:,3);[0,0,0]],'PlotStyle','compact')
%% discrimination between different pattern
datapool.R_betweenclust_pool = cell(5,1);
for j = 1:5
    datapool.R_betweenclust_pool{j} = [];
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{j})
            continue
        end
        for n = 1:length(data_decode(i).dataset{j})-1
            for m = n+1:length(data_decode(i).dataset{j})
            datapool.R_betweenclust_pool{j} = [datapool.R_betweenclust_pool{j};[data_decode(i).linear_regression_R{j}(n,m),data_decode(i).linear_regression_R_shuffle_clust{j}(n,m)]];
            end
        end
    end
end

datapool.error_betweenclust_pool = cell(5,1);
for j = 1:5
    datapool.error_betweenclust_pool{j} = [];
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{j})
            continue
        end
        for n = 1:length(data_decode(i).dataset{j})-1
            for m = n+1:length(data_decode(i).dataset{j})
            datapool.error_betweenclust_pool{j} = [datapool.error_betweenclust_pool{j};[data_decode(i).linear_regression_error{j}(n,m),data_decode(i).linear_regression_error_shuffle_clust{j}(n,m)]];
            end
        end
    end
end
idx_rmv = cell(1, 5);
for j = 1:5
    idx_rmv{j} = [idx_rmv{j},isnan(datapool.R_betweenclust_pool{j}(:,1))];
    idx_rmv{j} = [idx_rmv{j},isnan(datapool.R_betweenclust_pool{j}(:,2))];
    idx_rmv{j} = [idx_rmv{j},find(datapool.R_betweenclust_pool{j}(:,1)==0)];
    idx_rmv{j} = [idx_rmv{j},find(datapool.R_betweenclust_pool{j}(:,2)==0)];
    idx_rmv{j} = unique(idx_rmv{j});
    datapool.error_betweenclust_pool{j}(idx_rmv{j},:) = [];
    datapool.R_betweenclust_pool{j}(idx_rmv{j},:) = [];
end
%%
datapool.R_clustandspon_pool = cell(5,1);
datapool.idx_2strong = cell(5,1);
datapool.idx_1strong = cell(5,1);
datapool.idx_weak = cell(5,1);
datapool.idx_4clust = cell(5,1);
datapool.idx_3clust = cell(5,1);
datapool.idx_2clust = cell(5,1);
datapool.idx_1clust = cell(5,1);

for j = 1:5
    datapool.R_clustandspon_pool{j} = [];
    datapool.idx_2strong{j} = [];
    datapool.idx_1strong{j} = [];
    datapool.idx_weak{j} = [];
    datapool.idx_4clust{j} = [];
    datapool.idx_3clust{j} = [];
    datapool.idx_2clust{j} = [];
    datapool.idx_1clust{j} = [];
    count = 1;
    for i = 1:length(datapool.samecell_idx)
        if isempty(data_decode(i).correlation_matrix{j})
            continue
        end
        for n = 1:length(data_decode(i).dataset{j})
            if (isnan(data_decode(i).linear_regression_R_stim_spon{n,j}(2,1)))||(data_decode(i).linear_regression_R_stim_spon{n,j}(2,1)==0)
                continue
            end
            if (data_decode(i).strong_num(data_decode(i).dataset{j}(n))>=2)
                datapool.idx_2strong{j} = [datapool.idx_2strong{j}, count];
            elseif (data_decode(i).strong_num(data_decode(i).dataset{j}(n))==1)
                datapool.idx_1strong{j} = [datapool.idx_1strong{j}, count];
            else
                datapool.idx_weak{j} = [datapool.idx_weak{j}, count];
            end
            clust = data_decode(i).clust{data_decode(i).dataset{j}(n)};
            clust(find(clust==999)) = [];
            if length(clust)>=4
                datapool.idx_4clust{j} = [datapool.idx_4clust{j}, count];
            elseif length(clust)>=3
                datapool.idx_3clust{j} = [datapool.idx_3clust{j}, count];
            elseif length(clust)==2
                datapool.idx_2clust{j} = [datapool.idx_2clust{j}, count];
            else
                datapool.idx_1clust{j} = [datapool.idx_1clust{j}, count];
            end
%             a = data_decode(i).linear_regression_R_stim_spon{n,j}(2:6,1);
%             a(isnan(a)) = [];
%             datapool.R_clustandspon_pool{j}(count,1) = median(a);
%             [~,idx_shuffle] =min(abs(a-median(a)));
%             a = [data_decode(i).linear_regression_R_stim_spon{n,j}(2,3:6),data_decode(i).linear_regression_R_stim_spon{n,j}(3,4:6),data_decode(i).linear_regression_R_stim_spon{n,j}(4,5:6)];
%             a(isnan(a)) = [];
%             datapool.R_clustandspon_pool{j}(count,2) = a(idx_shuffle);
%             a = data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(2:6,1);
%             a(isnan(a)) = [];
%             datapool.R_clustandspon_pool{j}(count,3) = median(a);
%             [~,idx_shuffle] = min(abs(a-median(a)));
%             a = [data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(2,3:6),data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(3,4:6),data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(4,5:6)];
%             a(isnan(a)) = [];
%             datapool.R_clustandspon_pool{j}(count,4) = a(idx_shuffle);
            datapool.R_clustandspon_pool{j}(count,1) = data_decode(i).linear_regression_R_stim_spon{n,j}(2,1);
            datapool.R_clustandspon_pool{j}(count,3)  = data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(2,1);
            datapool.R_clustandspon_pool{j}(count,2)  = data_decode(i).linear_regression_R_stim_spon{n,j}(2,4);
            datapool.R_clustandspon_pool{j}(count,4)  = data_decode(i).linear_regression_R_stim_spon_shuffle{n,j}(2,4);
            count = count +1;
            
        end
    end
end
%%
g1 = repmat({'a'},size(datapool.idx_2strong{2},2),1);
g2 = repmat({'b'},size(datapool.idx_1strong{2},2),1);
g3 = repmat({'c'},size(datapool.idx_weak{2},2),1);
g4 = repmat({'d'},size(datapool.idx_4clust{2},2),1);
g5 = repmat({'e'},size(datapool.idx_3clust{2},2),1);
g6 = repmat({'f'},size(datapool.idx_2clust{2},2),1);
g7 = repmat({'g'},size(datapool.idx_1clust{2},2),1);
g11 = repmat({'a1'},size(datapool.idx_2strong{2},2),1);
g21 = repmat({'b1'},size(datapool.idx_1strong{2},2),1);
g31 = repmat({'c1'},size(datapool.idx_weak{2},2),1);
g41 = repmat({'d1'},size(datapool.idx_4clust{2},2),1);
g51 = repmat({'e1'},size(datapool.idx_3clust{2},2),1);
g61 = repmat({'f1'},size(datapool.idx_2clust{2},2),1);
g71 = repmat({'g1'},size(datapool.idx_1clust{2},2),1);
figure
y1 = [datapool.R_clustandspon_pool{2}(datapool.idx_2strong{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_2strong{2}, 2);datapool.R_clustandspon_pool{2}(datapool.idx_1strong{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_1strong{2}, 2);datapool.R_clustandspon_pool{2}(datapool.idx_weak{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_weak{2}, 2)];
boxplot(y1, [g1;g11;g2;g21;g3;g31],'Color', [datapool.colorpool(:,:,1);[0,0,0];datapool.colorpool(:,:,2);[0,0,0];datapool.colorpool(:,:,3);[0,0,0]],'PlotStyle','compact')
figure
y2 = [datapool.R_clustandspon_pool{2}(datapool.idx_4clust{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_4clust{2}, 2);datapool.R_clustandspon_pool{2}(datapool.idx_3clust{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_3clust{2}, 2);datapool.R_clustandspon_pool{2}(datapool.idx_2clust{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_2clust{2}, 2);datapool.R_clustandspon_pool{2}(datapool.idx_1clust{2}, 1);datapool.R_clustandspon_pool{2}(datapool.idx_1clust{2}, 2)];
boxplot(y2, [g4;g41;g5;g51;g6;g61;g7;g71],'Color', [datapool.colorpool(:,:,5);[0,0,0];datapool.colorpool(:,:,6);[0,0,0];datapool.colorpool(:,:,7);[0,0,0];datapool.colorpool(:,:,8);[0,0,0]],'PlotStyle','compact')
%%
a = [];
a(1) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_2strong{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_2strong{2}, 2))))/length(datapool.idx_2strong{2});
a(2) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_1strong{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_1strong{2}, 2))))/length(datapool.idx_1strong{2});
a(3) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_weak{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_weak{2}, 2))))/length(datapool.idx_weak{2});
a(4) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_4clust{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_4clust{2}, 2))))/length(datapool.idx_4clust{2});
a(5) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_3clust{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_3clust{2}, 2))))/length(datapool.idx_3clust{2});
a(6) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_2clust{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_2clust{2}, 2))))/length(datapool.idx_2clust{2});
a(7) = length(find((datapool.R_clustandspon_pool{2}(datapool.idx_1clust{2}, 1))>(datapool.R_clustandspon_pool{2}(datapool.idx_1clust{2}, 2))))/length(datapool.idx_1clust{2});
%%
g1 = repmat({1},size(datapool.R_pool{2},1),1);
g11 = repmat({2},size(datapool.R_pool{2},1),1);
g2 = repmat({3},size(datapool_ap5.R_pool{2},1),1);
g21 = repmat({4},size(datapool_ap5.R_pool{2},1),1);
g3 = repmat({5},size(datapool.R_pool{2},1),1);
g31 = repmat({6},size(datapool.R_pool{2},1),1);
g4 = repmat({7},size(datapool_ap5.R_pool{2},1),1);
g41 = repmat({8},size(datapool_ap5.R_pool{2},1),1);

y = [datapool.R_pool{2}(:,1);datapool.R_pool{2}(:,2);datapool_ap5.R_pool{2}(:,1);datapool_ap5.R_pool{2}(:,2);datapool.R_pool{4}(:,1);datapool.R_pool{4}(:,2);datapool_ap5.R_pool{4}(:,1);datapool_ap5.R_pool{4}(:,2)];
figure
boxplot(y, [g1;g11;g2;g21;g3;g31;g4;g41], 'Color',['r','k','b','k','r','k','b','k'], 'PlotStyle','compact')
%%
g1 = repmat({1},size(datapool.R_betweenclust_pool{2},1),1);
g11 = repmat({2},size(datapool.R_betweenclust_pool{2},1),1);
g2 = repmat({3},size(datapool_ap5.R_betweenclust_pool{2},1),1);
g21 = repmat({4},size(datapool_ap5.R_betweenclust_pool{2},1),1);
g3 = repmat({5},size(datapool.R_betweenclust_pool{4},1),1);
g31 = repmat({6},size(datapool.R_betweenclust_pool{4},1),1);
g4 = repmat({7},size(datapool_ap5.R_betweenclust_pool{4},1),1);
g41 = repmat({8},size(datapool_ap5.R_betweenclust_pool{4},1),1);

y = [datapool.R_betweenclust_pool{2}(:,1);datapool.R_betweenclust_pool{2}(:,2);datapool_ap5.R_betweenclust_pool{2}(:,1);datapool_ap5.R_betweenclust_pool{2}(:,2);datapool.R_betweenclust_pool{4}(:,1);datapool.R_betweenclust_pool{4}(:,2);datapool_ap5.R_betweenclust_pool{4}(:,1);datapool_ap5.R_betweenclust_pool{4}(:,2)];
figure
boxplot(y, [g1;g11;g2;g21;g3;g31;g4;g41], 'Color',['r','k','b','k','r','k','b','k'], 'PlotStyle','compact')
%%
g1 = repmat({1},size(datapool.corr_pool_byclust{2},1),1);
g11 = repmat({2},size(datapool.corr_pool_byclust{2},1),1);
g2 = repmat({3},size(datapool_ap5.corr_pool_byclust{2},1),1);
g21 = repmat({4},size(datapool_ap5.corr_pool_byclust{2},1),1);
g3 = repmat({5},size(datapool.corr_pool_byclust{4},1),1);
g31 = repmat({6},size(datapool.corr_pool_byclust{4},1),1);
g4 = repmat({7},size(datapool_ap5.corr_pool_byclust{4},1),1);
g41 = repmat({8},size(datapool_ap5.corr_pool_byclust{4},1),1);

y = [datapool.corr_pool_byclust{2}(:,1);datapool.corr_pool_byclust{2}(:,2);datapool_ap5.corr_pool_byclust{2}(:,1);datapool_ap5.corr_pool_byclust{2}(:,2);datapool.corr_pool_byclust{4}(:,1);datapool.corr_pool_byclust{4}(:,2);datapool_ap5.corr_pool_byclust{4}(:,1);datapool_ap5.corr_pool_byclust{4}(:,2)];
figure
boxplot(y, [g1;g11;g2;g21;g3;g31;g4;g41], 'Color',['r','k','b','k','r','k','b','k'], 'PlotStyle','compact')