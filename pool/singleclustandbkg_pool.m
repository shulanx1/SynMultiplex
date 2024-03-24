% compute width etc. for single cluster EPSP
datarange = [3,22];%1:length(data);

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

% compute spike time if not already computed
spike_thres =0.5;
thr = -10;
window_length = 0.1;
dt = 5e-5;

for i = datarange
    if ~isempty(data(i).stim_spike_time)
        continue
    end
    T = (0*dt:dt:(size(data(i).EPSP_highcond,1)-1)*dt);
    stim_start_new = 500;
    for j = 1:size(data(i).EPSP_highcond,2)
        dv_over_dt = diff(data(i).EPSP_highcond(:, j));
        idx_1 = [];
        for n1 = 2:size(data(i).EPSP_highcond,1)
            if data(i).EPSP_highcond(n1, j)> thr && data(i).EPSP_highcond(n1-1, j)<= thr
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
        data(i).stim_spike_time{j} = T(c(idx_2))-T(stim_start_new);
        data(i).stim_spike_time{j} = sort(data(i).stim_spike_time{j});
        data(i).stim_spike_time{j}(find((data(i).stim_spike_time{j}<0) |(data(i).stim_spike_time{j}>window_length))) = [];
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
    if (isempty(data(i).jitter))||(isnan(data(i).jitter))
        a = [];
        fired_trial = 0;
        burst_trial = 0;
        for j = 1:size(data(i).stim_spike_time,2)
            if ~isempty(data(i).stim_spike_time{j})
                a = [a, min(data(i).stim_spike_time{j})];
                fired_trial = fired_trial + 1;
                if length(data(i).stim_spike_time{j})>=2
                    burst_trial = burst_trial + 1;
                end
            end
        end
        data(i).p_spike = fired_trial/size(data(i).stim_spike_time,2);
        data(i).p_burst = burst_trial/size(data(i).stim_spike_time,2);
        if size(data(i).EPSP_highcond,2)>1
            data(i).jitter = std(a*1e3);
        end
    end
    data(i).mean_onset = 0;
    data(i).onset = 0;
    data(i).median_onset = 0;
    data(i).AP_num = 0;

    data(i).isi = [];
    if isempty(cell2mat(data(i).stim_spike_time))
        continue
    else
        a = [];
        for n = 1:size(data(i).stim_spike_time,2)
            if ~isempty(data(i).stim_spike_time{n})
                a = [a,min(data(i).stim_spike_time{n})];
            end
        end
        data(i).mean_onset = mean(a)*1e3;
        data(i).median_onset = median(a)*1e3;
        data(i).onset = min(cell2mat(data(i).stim_spike_time))*1e3;
        data(i).AP_num = length(cell2mat(data(i).stim_spike_time))/10;
        for n = 1:size(data(i).stim_spike_time,2)
            if length(data(i).stim_spike_time{n})<=1
                continue
            else
                data(i).isi = [data(i).isi, diff(data(i).stim_spike_time{n})];
            end
        end
    end

end
for i = datarange
    
    if isempty(data(i).stim_spike_time)
        continue
    end
    if isempty(cell2mat(data(i).stim_spike_time))
        continue
    else
        p_burst = 0;
        for n = 1:size(data(i).stim_spike_time,2)
            if length(data(i).stim_spike_time{n})<=1
                continue
            else                    
                if min(diff(data(i).stim_spike_time{n})*1e3)<= 20
                    p_burst = p_burst + 1;
                end
            end
        end
        data(i).p_burst_true = p_burst/size(data(i).stim_spike_time,2);
    end
end
%% pool
datapool.idx_strong = [];
datapool.idx_weak = [];
for i = 1:length(data)
    if data(i).maxdv_single(2)>2
        datapool.idx_strong = [datapool.idx_strong,i];
    else
        datapool.idx_weak = [datapool.idx_weak,i];
    end
end
datapool.jitter_strong = [data(datapool.idx_strong).jitter];
datapool.jitter_strong(datapool.jitter_strong==0) = [];
datapool.jitter_weak = [data(datapool.idx_weak).jitter];
datapool.jitter_weak(datapool.jitter_weak==0) = [];

datapool.onset_strong = [data(datapool.idx_strong).onset];
datapool.onset_strong(datapool.onset_strong==0) = [];
datapool.onset_weak = [data(datapool.idx_weak).onset];
datapool.onset_weak(datapool.onset_weak==0) = [];

datapool.idx_burst = [];
datapool.idx_single = [];
for i = 1:length(data)
    if (data(i).p_burst>0)&&(data(i).p_spike>0)
        datapool.idx_burst = [datapool.idx_burst,i];
    elseif data(i).p_spike>0
        datapool.idx_single = [datapool.idx_single,i];
    end
end
datapool.area_burst = [];
datapool.width_burst = [];
for i = datapool.idx_burst
    datapool.area_burst = [datapool.area_burst, data(i).area_single(2)];
    datapool.width_burst = [datapool.width_burst, data(i).halfwidth_single(2)];
end
datapool.area_single = [];
datapool.width_single = [];
for i = datapool.idx_single
    datapool.area_single = [datapool.area_single, data(i).area_single(2)];
    datapool.width_single = [datapool.width_single, data(i).halfwidth_single(2)];
end