
idx_single =[7];
cluster = {[1:4]};



path = 'E:\data\uncaging\multi cluster cooperation\112222\cell4';
preflix = '221123_001';

% single clust
clear PC shutter 
filename_single = sprintf('%s.singlecluster.%d.wcp',preflix, idx_single);
out=import_wcp(fullfile(path, filename_single),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP_single = [];
figure
for i = 1:n_recording
    Vm_single(:,i) = out.S{3}(:,i);
    Im_single(:,i)  = out.S{4}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>50);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>1000)+1]);
    for j = 1:length(stim_start{i})
        if shutter(stim_start{i}(j), i)< 4
            continue
        else
            stim_start_new = 500;
            EPSP_single(:,j) = Vm_single(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i);%-Vm_SLM(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp_single(j),max_amp_single(j)] = max(abs(EPSP_single(stim_start_new:stim_start_new+1200,j)-EPSP_single(stim_start_new,j)));
            max_amp_single(j) = max_amp_single(j)+stim_start_new;
            [half,I] = sort(abs(EPSP_single(:,j)-amp_single(j)/2));
            b = find(I<max_amp_single(j) & I>stim_start_new);
            a = find(I>max_amp_single(j));
            halfwidth_single(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_single(j) = (max_amp_single(j)-stim_start_new)*dt*1000; %in ms
            maxdv_single(j) = max(diff(movmean(EPSP_single(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_single(:,j), 10)));
            subplot(1,5,j)
            plot(t, EPSP_single(:,j), 'Color', [0.5,0.5,0.5])
            hold on
            scatter(t([I(a(1)),I(b(1)),max_amp_single(j)]), EPSP_single([I(a(1)),I(b(1)),max_amp_single(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
            scatter(t(c), EPSP_single(c, j))
            xlabel('t (ms)')
            ylabel('Vm (mV)')
            xlim([min(t), max(t)])
            box off
        end
    end
end
Vm_single(stim_start{1}(1), 1)
save(fullfile(path, sprintf('singleclust_%d.mat', idx_single)),'EPSP_single', 'Vm_single', 'amp_single', 'risetime_single','halfwidth_single', 'maxdv_single','t','idx_single','dt')
%% multiclust_gradual_sum
idx_multi =[2];



clear PC shutter
filename_multi = sprintf('%s.multicluster_sum.%d.wcp',preflix, idx_multi);
out=import_wcp(fullfile(path, filename_multi),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;
clear PC shutter
EPSP_multi = [];
figure
for i = 1:n_recording
    Vm_multi(:,i) = out.S{3}(:,i);
    Im_multi(:,i)  = out.S{4}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>50);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>1000)+1]);
    for j = 1:length(stim_start{i})
        if shutter(stim_start{i}(j), i)< 4
            continue
        else
            stim_start_new = 500;
            EPSP_multi(:,j) = Vm_multi(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i);%-Vm_SLM(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp_multi(j),max_amp_multi(j)] = max(abs(EPSP_multi(stim_start_new:stim_start_new+1200,j)-EPSP_multi(stim_start_new,j)));
            max_amp_multi(j) = max_amp_multi(j)+stim_start_new;
            [half,I] = sort(abs(EPSP_multi(:,j)-amp_multi(j)/2));
            b = find(I<max_amp_multi(j) & I>stim_start_new);
            a = find(I>max_amp_multi(j));
            halfwidth_multi(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_multi(j) = (max_amp_multi(j)-stim_start_new)*dt*1000; %in ms
            maxdv_multi(j) = max(diff(movmean(EPSP_multi(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_multi(:,j), 10)));
            subplot(1,5,j)
            plot(t, EPSP_multi(:,j), 'Color', [0.5,0.5,0.5])
            hold on
            scatter(t([I(a(1)),I(b(1)),max_amp_multi(j)]), EPSP_multi([I(a(1)),I(b(1)),max_amp_multi(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
            scatter(t(c), EPSP_multi(c, j))
            xlabel('t (ms)')
            ylabel('Vm (mV)')
            xlim([min(t), max(t)])
            box off
        end
    end
end
Vm_multi(stim_start{1}(1), 1)
save(fullfile(path, sprintf('multiclust_sum_%d.mat', idx_multi)),'EPSP_multi', 'Vm_multi', 'amp_multi', 'risetime_multi','halfwidth_multi', 'maxdv_multi','t','idx_multi','dt')

%% highcond
close all
idx_highcond = 4;


filename = sprintf('%s.multicluster_highcond.%d.wcp', preflix, idx_highcond);
out=import_wcp(fullfile(path, filename),'debug');

n_channel = out.channel_no;
n_recording = length(out.rec_index);
clear Vm IM PC shutter isi mean_FR
stim_start = cell(1,n_recording);
EPSP = cell(1, n_recording);
V_reststate = [];
jitter = [];
p_spike = [];

T = out.T;
dt = T(2) - T(1);
t = (-500*dt:dt:5500*dt)*1e3;
for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>100);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_eachclust = find(diff(PC(:,i))>100);
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>1000)+1]);
    for j = 1:length(stim_start{i})
        if shutter(stim_start{i}, i)< 4
            continue
        else
        EPSP{i} = [EPSP{i}, Vm(stim_start{i}(j)-500:stim_start{i}(j)+5500,i)];
        RMP(i) = Vm(stim_start{i}(j), i);
        end
    end
    V_reststate = [V_reststate;Vm(1:floor(0.3/dt) ,i)-RMP(i)];
end


% figure
% plot(t, EPSP, 'Color', [0.5,0.5,0.5], 'LineWidth', 0.5)
% hold on
% plot(t,mean_EPSP, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2)
% hold on
% scatter(t([I(a(1)),I(b(1)),max_amp]), mean_EPSP([I(a(1)),I(b(1)),max_amp],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
% hold on
% scatter(t(c), mean_EPSP(c))
% xlabel('t (ms)')
% ylabel('Vm (mV)')
% xlim([min(t), max(t)])

% box off
figure
for i = 1:n_recording
    subplot(1,5,i)
    plot(t, EPSP{i},'k', 'LineWidth', 1)
    hold on
    y1 = ylim();
    yplot = [y1(1):0.1:y1(2)];
    for i = 1:round(length(    stim_eachclust)/10)
        plot(ones(1, length(yplot))*t(stim_eachclust(i)-stim_eachclust(1)+500),yplot, 'Color', [0.5,0.5,0.5])
    end
end
% axis tight
xlim([-20,200])
rms_vrest = rms(V_reststate);
meanRMP = mean(RMP);

spike_thres =0.6;
thr = 0;
stim_spike_time = cell(n_recording, length(stim_start{1}));
window_length = 0.1; 
for i = 1:n_recording
    dv_over_dt = diff(Vm(:,i));
%     [~,b] = findpeaks(Vm(:,i));
%     idx_1 = find(Vm(b,i)>=0);
%     idx_1(find(diff(idx_1)<=1)+1)=[];
    idx_1 = [];
    for n = 2:size(Vm,1)
        if Vm(n, i)> thr && Vm(n-1, i)<= thr
            idx_1 = [idx_1, n]; % spike time in s
        end
    end
    n_spike = length(idx_1);
    ddv_over_dt = diff(dv_over_dt);
    [~,c] = findpeaks(dv_over_dt);
    idx = find(dv_over_dt(c)>=spike_thres);
    idx(find(diff(idx)<=2)+1) = [];

    idx_2 = [];
    for n = 1:length(idx_1)
        a = abs(c(idx)-idx_1(n));
        I = find(a==min(a));
        I = I(1);
        idx_2(n) = idx(I);
        idx(I) = [];
    end
    spike_point{i} = c(idx_2);
    spike_time{i} = T(c(idx_2));
    spike_time{i} = sort(spike_time{i});
    FR(i) = length(spike_time{i}>=10)/5;
    stim_start_time{i} = T(stim_start{i});
for j = 1:length(stim_start{i})
    fired_trial(i, j) = ~isempty(find((spike_time{i}>=stim_start_time{i}(j))&(spike_time{i}<=stim_start_time{i}(j)+window_length)));
    burst_trial(i, j) = length(find((spike_time{i}>=stim_start_time{i}(j))&(spike_time{i}<=stim_start_time{i}(j)+window_length)))>=2;
    if ~isempty(find((spike_time{i}>=stim_start_time{i}(j))&(spike_time{i}<=stim_start_time{i}(j)+window_length)))
        a = find((spike_time{i}>=stim_start_time{i}(j))&(spike_time{i}<=stim_start_time{i}(j)+window_length));
        stim_spike_time{i, j} = spike_time{i}(a)-stim_start_time{i}(j);
    end
end
    isi{i} = diff(spike_time{i});
end

isi = cell2mat(isi);
mean_FR = mean(FR);
P_spike= [];
P_burst = [];
jitter = [];
onset = [];

for i = 1:n_recording
    P_spike(i) = length(find(fired_trial(i,:))~=0)/10;
    P_burst(i) = length(find(burst_trial(i,:))~=0)/10;%(P_spike(i)*10);
    a = [];
    for j = 1:length(stim_start_time{i})
        if ~isempty(stim_spike_time{i, j})
            a = [a, min(stim_spike_time{i, j})];
        end
    end
    if ~isempty(a)
        jitter(i) = std(a*1000);%(max(stim_spike_time)-min(stim_spike_time))*1000
        onset(i) = min(a)*1000;
    end
end

FF = std(FR)^2/mean(FR);
CV = std(isi)/mean(isi);

% for i =1:length(spike_time)
%     spike_time1{i} = 1000*(spike_time{i}-0.5);
% end
% figure
% %histogram((spike_time1-0.5)*1000, [-20:5:200]);
% [xPoints, yPoints] = plotSpikeRaster(spike_time1, 'PlotType', 'vertline');
% hold on
% y1 = ylim();
% yplot = [y1(1):0.1:y1(2)];
% for i = 1:length(stim_eachclust)
%     plot(ones(1, length(yplot))*t(stim_eachclust(i)-stim_eachclust(1)+500),yplot, 'Color', [0.5,0.5,0.5])
% end
% xlim([-10,100])

figure
for i = 1:n_recording
    if isempty(spike_time{i})
        continue
    end
    subplot(1, 5, i)
    spike_time1 = cell(1, length(stim_start_time{i}));
    for j = 1:length(stim_start_time{i})
        spike_time1{j} = (spike_time{i} - stim_start_time{i}(j))*1000;
    end
    [xPoints, yPoints] = plotSpikeRaster(spike_time1, 'PlotType', 'vertline');
    hold on
    y1 = ylim();
    yplot = [y1(1):0.1:y1(2)];
    for i = 1:length(stim_eachclust)/10
        plot(ones(1, length(yplot))*t(stim_eachclust(i)-stim_eachclust(1)+500),yplot, 'Color', [0.5,0.5,0.5])
    end
    xlim([-20,200])
end
out_file = fullfile(path, sprintf('multiclust_highcond_%d.mat', idx_highcond));
save(out_file, 'RMP','Vm', 'EPSP', 'stim_start', 't', 'T', 'stim_start_time', 'stim_eachclust', 'spike_time', 'stim_spike_time', 'fired_trial', 'burst_trial', 'P_spike', 'P_burst', 'jitter', 'onset');



%%
% linear_regression_R_PCA = cell(5, length(samecell_idx));
% linear_regression_error_PCA = cell(5, length(samecell_idx));
linear_regression_R_shuffle_clust = cell(5, length(samecell_idx));
linear_regression_error_shuffle_clust = cell(5, length(samecell_idx));
for i = 1:length(samecell_idx)
    for j = 1:5        
         A = cell2mat(stim_spike_time_matrix_percell{j, i}');
%         [linear_regression_R_PCA{j,i}, linear_regression_error_PCA{j,i}] = multicluster_linear_regression(A, n_trials);
    A_shuffle = A(randperm(size(A, 1)), :);
    [linear_regression_R_shuffle_clust{j,i}, linear_regression_error_shuffle_clust{j,i}] = multicluster_linear_regression(A_shuffle, n_trials);
    for n = 1:9
        A_shuffle = A(randperm(size(A, 1)), :);
        [R_shuffle_clust, error_shuffle_clust] = multicluster_linear_regression(A_shuffle, n_trials);
        linear_regression_R_shuffle_clust{j,i} = linear_regression_R_shuffle_clust{j,i} + R_shuffle_clust;
        linear_regression_error_shuffle_clust{j,i} = linear_regression_error_shuffle_clust{j,i} + error_shuffle_clust;
    end
    linear_regression_R_shuffle_clust{j,i} = linear_regression_R_shuffle_clust{j,i}/10;
    linear_regression_error_shuffle_clust{j,i} = linear_regression_error_shuffle_clust{j,i}/10;
    end
end
%%
R_pool = cell(5,1);
for n = 1:5
        R_pool{n} = zeros(length(samecell_idx), 2);
    for i = 1:length(samecell_idx)
        R_pool{n}(i,1) = max(max(linear_regression_R_PCA{n,i}));
        R_pool{n}(i,2) = max(max(linear_regression_R_shuffle_clust{n,i}));
    end
end
error_pool = cell(5,1);
for n = 1:5
        error_pool{n} = zeros(length(samecell_idx), 2);
    for i = 1:length(samecell_idx)
        error_pool{n}(i,1) = max(max(linear_regression_error_PCA{n,i}));
        error_pool{n}(i,2) = max(max(linear_regression_error_shuffle_clust{n,i}));
    end
end