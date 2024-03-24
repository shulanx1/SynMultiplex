% datarange = 1:length(data);
datarange = 137:139;
dt = 5e-5;
for i =datarange
    if isempty(data(i).EPSP)
        continue
    end
    data(i).amp_single = [];
    data(i).halfwidth_single = [];
    data(i).maxdv_single = [];
    data(i).area_single = [];
    data(i).risetime_single = [];
    if isempty(data(i).EPSP)
        continue
    end
    for j = 1:size(data(i).EPSP,2)
        stim_start_new = 500;
        [data(i).amp_single(j),max_amp] = max(abs(data(i).EPSP(stim_start_new:stim_start_new+2000, j)-data(i).EPSP(stim_start_new, j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs((data(i).EPSP(:,j)-data(i).EPSP(stim_start_new, j))-data(i).amp_single(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        data(i).halfwidth_single(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
        data(i).risetime_single(j) = (max_amp-stim_start_new)*dt*1000; %in ms
        data(i).maxdv_single(j) = max(diff(movmean(data(i).EPSP(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
        EPSP_mean = data(i).EPSP(:,j)-data(i).EPSP(stim_start_new, j);
        data(i).area_single(j) = sum(EPSP_mean(501:2501))*dt;
    end
end

for i =datarange
    if isempty(data(i).EPSP_singlespine)
        continue
    end

    data(i).amp_spine = [];
    data(i).halfwidth_spine = [];
    data(i).maxdv_spine = [];
    data(i).area_spine = [];
    data(i).risetime_spine = [];
    if isempty(data(i).EPSP_singlespine)
        continue
    end
    for j = 1:size(data(i).EPSP_singlespine,2)
        stim_start_new = 500;
        [data(i).amp_spine(j),max_amp] = max(abs(data(i).EPSP_singlespine(stim_start_new:stim_start_new+2000, j)-data(i).EPSP_singlespine(stim_start_new, j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs((data(i).EPSP_singlespine(:,j)-data(i).EPSP_singlespine(stim_start_new, j))-data(i).amp_spine(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        data(i).halfwidth_spine(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
        data(i).risetime_spine(j) = (max_amp-stim_start_new)*dt*1000; %in ms
        data(i).maxdv_spine(j) = max(diff(movmean(data(i).EPSP_singlespine(stim_start_new:stim_start_new+1000, j), 10)))/dt/1000; %in V/s
        EPSP_mean = data(i).EPSP_singlespine(:,j)-data(i).EPSP_singlespine(stim_start_new, j);
        data(i).area_spine(j) = sum(EPSP_mean(501:2501))*dt;
    end

end

%%
datapool.strong_idx = [];
datapool.weak_idx = [];
for i = 1:length(data)
    if (max(data(i).maxdv_single)>2)&(max(data(i).maxdv_single)<20)
        datapool.strong_idx =[ datapool.strong_idx,i];
    elseif (max(data(i).maxdv_single)<20)
        datapool.weak_idx =[ datapool.weak_idx,i];
    end
end

datapool.weak_amp = [];
datapool.strong_amp = [];
count = 1;
for i = datapool.strong_idx
    if max(data(i).amp_single)<30
        datapool.strong_amp(count,:) = data(i).amp_single(end-6:end);
        count = count +1;
    end
end

count = 1;
for i = datapool.weak_idx
    datapool.weak_amp(count,:) = data(i).amp_single(end-6:end);
    count = count +1;
end

datapool.strong_amp_1spine = [];
for i = datapool.strong_idx
    if size(data(i).EPSP,2)==8
        datapool.strong_amp_1spine = [datapool.strong_amp_1spine,min(data(i).amp_single(1),data(i).amp_spine(1),data(i).amp_single(2))];
    elseif ~isempty(data(i).EPSP_singlespine)
        datapool.strong_amp_1spine = [datapool.strong_amp_1spine,min(data(i).amp_spine(1), data(i).amp_spine(2))];
    end
end

datapool.weak_amp_1spine = [];
for i = datapool.weak_idx
    if ~isempty(data(i).EPSP_singlespine)
        datapool.weak_amp_1spine = [datapool.weak_amp_1spine,min(data(i).amp_spine(1), data(i).amp_spine(2))];
    end
end
%%
datapool.weak_large_amp = [];
count = 1;
for i = 1:size(datapool.weak_amp,1)
    if (max(datapool.weak_amp(i,:)>8))&(datapool.weak_amp(i,4)<6)
        datapool.weak_large_amp(count,:) = datapool.weak_amp(i,:);
        count = count +1;
    end
end
%%
datapool.weak_norm = [];
datapool.strong_norm = [];
for i = 1:size(datapool.weak_large_amp,1)
    datapool.weak_norm(i,:) = datapool.weak_large_amp(i,:)/max(datapool.weak_large_amp(i,:));
end
for i = 1:size(datapool.strong_amp,1)
    datapool.strong_norm(i,:) = datapool.strong_amp(i,:)/max(datapool.strong_amp(i,:));
end