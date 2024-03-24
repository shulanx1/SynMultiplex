% close all
clear all


idx_SLM = 2;
branch_no = 2;


path = 'E:\data\uncaging\multi cluster cooperation\010322\cell3';
preflix = '220103_001';


clear PC shutter
filename_SLM = sprintf('%s.summation_rep.%d.wcp',preflix, idx_SLM);
out=import_wcp(fullfile(path, filename_SLM),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);
%  n_recording =  4;

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    Im(:,i) = out.S{4}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>50);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>50)+1]);
    stim_start_time{i} = T(stim_start{i});
end

for j = 1:length(stim_start{1})
    Vm_seg{j} = [];
    for i = 1:n_recording
        Vm_seg{j} = [Vm_seg{j},Vm(stim_start{i}(j)-500:stim_start{i}(j)+4500,i)];
    end
end
figure
plot(Vm)


spike_thres = 0.35;
thr = 0;
stim_spike_time = [];

% for i = 1:n_recording
%     dv_over_dt = diff(Vm(:,i));
%     [~,b] = findpeaks(Vm(:,i));
%     idx_1 = find(Vm(b,i)>=0);
% %      idx_1(find(diff(idx_1)<=4)+1)=[];
%     n_spike = length(idx_1);
%     ddv_over_dt = diff(dv_over_dt);
%     [~,c] = findpeaks(ddv_over_dt);
%     idx = find(ddv_over_dt(c)>=spike_thres);
%     idx(find(diff(idx)<=2)+1) = [];
%     [~,I] = sort(ddv_over_dt(c(idx)));
%     idx = idx(I(end-n_spike+1:end));
%     spike_point{i} = c(idx);
%     spike_time{i} = T(c(idx));
% %     spike_hist(i, c(idx)) = 1;   
%     spike_time{i} = sort(spike_time{i});
% %     for k = 1:length(max_point)
% %         for m = max_point(k):-1:max_point(k)-30
% %             if ddv_over_dt(m)*ddv_over_dt(m-1)<=0
% %                 spike_point{i}(k) = m;
% %                 spike_time{i}(k) = T(m);
% %                 break
% %             end
% %         endspike_T
% %     end
%     FR(i) = length(find(spike_time{i}>=7))/2;
%     for j = 1:length(stim_start_time{i})
% %         spike_time_trial{i,j} = spike_time{i}(find(spike_time{i}>=stim_start_time{i}(j)&spike_time{i}<=stim_start_time{i}(j)+0.05));
% spike_time_trial{i,j} = spike_time{i}(find(spike_time{i}>=stim_start_time{i}(j)&spike_time{i}<=stim_start_time{i}(j)+0.1));
%     end
% end
for i = 1:n_recording
    dv_over_dt = diff(Vm(:,i));
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
    FR(i) = length(find(spike_time{i}>=17))/3;
end

T_down = 0:0.001:T(end);
T_down = T_down(1:end-1)+0.0005;
psth_down = hist(cell2mat(spike_time),T_down);
mean_FR = mean(FR);
FR_dynamic = adaptive_window_length(psth_down, T_down, min(max(3,round(mean(FR))),10))/10;
[~, a] = find(FR_dynamic>=max(mean_FR,0.01));
d = find(diff(a)~=1);
event_start = [a(1), a(d+1)];
event_end = [a(d), a(end)];
for i = 1:length(stim_start_time{1})
    stim_diff = abs(T_down(event_start)-stim_start_time{1}(i));
    if min(stim_diff)>=0.4
        FR_temp(i) = 0;
        FR_area(i) = 0;
        FR_ave(i) = 0;
        continue
    end    
    a = find(stim_diff==min(stim_diff));
    FR_temp(i) = max(FR_dynamic(event_start(a):event_end(a)));
    FR_area(i) = sum(FR_dynamic(event_start(a):event_end(a)))*0.001;
    FR_ave(i) = mean(FR_dynamic(event_start(a):event_end(a)));
    event_start_stim(i) = round(T_down(event_start(a)))/dt;
    event_end_stim(i) = round(T_down(event_end(a))/dt);
end
num = zeros(1,length(stim_start{1}));
precision = [];
spike_within_event = {};
for j = 1:length(event_start_stim)
    t_event = [];
    for i = 1:n_recording
        spike_within_event{i,j} = find(spike_point{i}>=event_start_stim(j) & spike_point{i}<=event_end_stim(j));
        if ~isempty(spike_within_event{i,j})
            num(j) = num(j)+1;
            t_event = [t_event; dt*spike_within_event{i,j}];
        end
    end
    precision(j) = std(t_event);
end
reliability = num/10;
p_spike = zeros(1,size(spike_time_trial,2));
for i = 1:size(spike_time_trial,2)
for j = 1:size(spike_time_trial,1)
if ~isempty(spike_time_trial{j,i})
p_spike(i) = p_spike(i)+1;
end
end
end
figure
plot(T_down, FR_dynamic)
p_spike = p_spike/10;
save(fullfile(path, sprintf('summation_rep_%d.mat', idx_SLM)),'Vm','Vm_seg','Im', 'T','idx_SLM','dt', 'stim_start', 'FR_ave', 'FR_area','FR_temp','stim_start_time', 'spike_time','reliability','precision', 'FR_dynamic', 'FR', 'spike_time_trial', 'mean_FR','p_spike')
%%
% spike_time1 = sort(cell2mat(spike_time));
% for i = 1:7
% spike_time_trial1{i} = sort(cell2mat(spike_time_trial(:,i)'));
% end
% spike_time_trial_total = cell2mat(spike_time_trial1);
% [~,idx]=intersect(spike_time1,spike_time_trial_total);
% spike_time1(idx) = [];
% mean_isi = mean(diff(spike_time1));
% mean_FR1 = 10/mean_isi;
% for i = 1:7
%     stim_isi(i) = mean(diff(spike_time_trial1{i}));
% end
% stim_FR1 = 10./stim_isi;
% spike_time1 = spike_time;
% for j = 1:size(spike_time_trial,1)
% for i = 1:size(spike_time_trial,2)
% if ~isempty(spike_time_trial{j,i})
%     [~,idx]=intersect(spike_time1{j},spike_time_trial{j,i});
%     a = diff(spike_time1{j});
%     if isempty(find(idx>1))
%         continue
%     end
%     isi_stim(j,i) = mean(a(idx-1));
%     spike_time1{j}(idx) = [];
%     clear a
% end
% isi(j) = mean(diff(spike_time1{j}));
% end
% end
% mean_isi = mean(isi);
% mean_FR1 = 1/mean_isi;
% for i = 1:7
% if isempty(find(isi_stim(:,i)~=0))
% stim_FR1(i) = 0;
% continue
% else
% idx1 = find(isi_stim(:,i)~=0);
% stim_FR1(i) = 1/mean(isi_stim(idx1,i));
% end
% end
for i = 1:size(spike_time_trial,1)
for j = 1:length(stim_start_time{i})
%         spike_time_trial{i,j} = spike_time{i}(find(spike_time{i}>=stim_start_time{i}(j)&spike_time{i}<=stim_start_time{i}(j)+0.05));
spike_time_trial1{i,j} = spike_time{i}(find(spike_time{i}>stim_start_time{i}(j)-0.1&spike_time{i}<=stim_start_time{i}(j)+0.9));
[~,idx]=intersect(spike_time{i},spike_time_trial1{i,j});
a = diff(spike_time{i});
if ~isempty(find(idx==1))
    idx(find(idx==1))=[];
end
isi1{i,j} = a(idx-1);
FR1(i, j) = length(spike_time_trial1{i,j});
end
spike_time_control{i} = spike_time{i}(find(spike_time{i}>7));
[~,idx]=intersect(spike_time{i},spike_time_control{i});
a = diff(spike_time{i});
isi1{i,8} = a(idx-1);
FR1(i, 8) = length(spike_time_control{i})/2;
end
for i = 1:8
%     isi_vector{i} = cell2mat(isi1(:,i)');
%     CV(i) = std(isi_vector{i})/mean(isi_vector{i});
    FF(i) = std(FR1(:,i))^2/mean(FR1(:,i));
end
FF([1,end-2,end])
%% psth
edges = [0:0.1:9];
figure
hist(cell2mat(spike_time), edges)
%% pool
AP_thre = [];
RMP_thre = [];
onset_delay = {};
spike_jitter = [];
p_burst = [];
spike_thres = 0.35;
T = (-500*dt:dt:4500*dt)*1e3;
RMP = {};

% stim_response1 = stim_response;
% stim_response = reshape(stim_response, 1,[]);
% for i = 1:length(stim_response)
%     stim_response{i} = stim_response{i} - stim_response{i}(500,:);
% end
for i = 1:length(stim_response)
spike_time = {};
first_spike_time = [];
    if isempty(find(stim_response{i})>=50)
        continue
    end
    n_spike_trial = [];
    burst_trial = [];
    first_spike_time = [];
    for j = 1:size(stim_response{i}, 2)
        dv_over_dt = diff(stim_response{i}(:,j)');
        [~,b] = findpeaks(stim_response{i}(:,j)');
        idx_1 = find(stim_response{i}(b,j)>=50);
        idx_2 = b(idx_1);
        n_spike = length(idx_1);
        idx_2(find(idx_2<=500))=[];
        n_spike_trial = [n_spike_trial, n_spike];
        if ~isempty(find(diff(idx_2)<=500))
            burst_trial = [burst_trial,j];
        end
        a = sort(idx_2);
        spike_time{j} = T(a); 
        first_spike_time = [first_spike_time, min(T(a))];
        clear a
    end
    

    onset_delay{i} = first_spike_time;
    a = first_spike_time;
    if isempty(a)
        continue
    end
%     spike_jitter(i) = std(a);
    p_burst(i) = length(burst_trial)/length(find(n_spike_trial>0));
end
%%
%  onset_delay(34:36) = [0,0,0];
% p_burst(34:36) = [0,0,0];
% p_burst = reshape(p_burst,9,[]);
% onset_delay = reshape(onset_delay,9,[]);
y1 = cell2mat(onset_delay(:,1)');
y2 = cell2mat(onset_delay(:,2)');
y3 = cell2mat(onset_delay(:,3)');
y4 = cell2mat(onset_delay(:,4)');
g1 = repmat({'a'},length(y1),1);
g2 = repmat({'b'},length(y2),1);
g3 = repmat({'c'},length(y3),1);
g4 = repmat({'d'},length(y4),1);
g = [g1;g2;g3;g4];
y = [y1';y2';y3';y4'];
len = 4;
pink = [1, 0, 0];
red = [1,0.7,0.7];
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
figure
boxplot(y,g, 'Color', colors_p,'PlotStyle','compact')
%%
%%
close all
% clear all


idx_SLM = 16;
branch_no = 1;


path = 'E:\data\uncaging\multi cluster cooperation\010322\cell3';
preflix = '220103_001';


clear PC shutter Vm
filename_SLM = sprintf('%s.summation_rep.%d.wcp',preflix, idx_SLM);
out=import_wcp(fullfile(path, filename_SLM),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);
% n_recording =  8;

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    Im(:,i) = out.S{4}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>100);
    if size(stim_start{i}, 1) == 0
        continue
    end
    idx1 = [1; find(diff(stim_start{i})>100)+1];
    stim_start{i} = stim_start{i}(idx1);
    stim_start_time{i} = T(stim_start{i});
end
stim_start_1 = cell2mat(stim_start');
stim_start_time_1 = T(stim_start_1);
n_stim  = 7;
n_trials = length(stim_start_1)/n_stim;

for j = 1:n_stim
    Vm_seg{j} = [];
    idx1 = [0:n_trials-1]*n_stim + j;
    for i = idx1
        stim_idx = mod(i, length(stim_start{1}));
        if stim_idx == 0
            stim_idx = length(stim_start{1});
        end
        recording_idx = ceil(i/length(stim_start{1}));
        Vm_seg{j} = [Vm_seg{j},Vm(stim_start{recording_idx}(stim_idx)-500:stim_start{recording_idx}(stim_idx)+4500,recording_idx)];
    end
end




spike_thres = 1;
thr = 0;
stim_spike_time = [];
% for i = 1:n_recording
%     dv_over_dt = diff(Vm(:,i));
%     [~,b] = findpeaks(Vm(:,i));
%     idx_1 = find(Vm(b,i)>=0);
%     n_spike = length(idx_1);
%     ddv_over_dt = diff(dv_over_dt);
%     [~,c] = findpeaks(ddv_over_dt);
%     idx = find(ddv_over_dt(c)>=spike_thres);
%     idx(find(diff(idx)<=2)+1) = [];
%     [~,I] = sort(ddv_over_dt(c(idx)));
%     idx = idx(I(end-n_spike+1:end));
%     spike_point{i} = c(idx);
%     spike_time{i} = T(c(idx));
%   
%     spike_time{i} = sort(spike_time{i});
%     FR(i) = length(find(spike_time{i}>=17))/3;
% end
for i = 1:n_recording
    dv_over_dt = diff(Vm(:,i));
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
    FR(i) = length(find(spike_time{i}>=17))/3;
end

window_length = 0.03;
spike_time_trial = cell( n_trials, n_stim);
for j = 1:n_stim
    idx1 = [0:n_trials-1]*n_stim + j;
    for i = idx1
        stim_idx = mod(i, length(stim_start{1}));
        if stim_idx == 0
            stim_idx = length(stim_start{1});
        end
        recording_idx = ceil(i/length(stim_start{1}));
        trial_idx = find(i==idx1);
        spike_time_trial{trial_idx,j} = spike_time{recording_idx}(find(spike_time{recording_idx}>=stim_start_time_1(i)&spike_time{recording_idx}<=stim_start_time_1(i)+window_length))- stim_start_time_1(i);
    end
end
mean_FR = mean(FR);
% T_down = 0:0.001:T(end);
% T_down = T_down(1:end-1)+0.0005;
% psth_down = hist(cell2mat(spike_time),T_down);
% mean_FR = mean(FR);
% FR_dynamic = adaptive_window_length(psth_down, T_down, min(max(3,round(mean(FR))),10))/10;
% [~, a] = find(FR_dynamic>=max(mean_FR,0.01));
% d = find(diff(a)~=1);
% event_start = [a(1), a(d+1)];
% event_end = [a(d), a(end)];
% for i = 1:length(stim_start_time{1})
%     stim_diff = abs(T_down(event_start)-stim_start_time{1}(i));
%     if min(stim_diff)>=0.4
%         FR_temp(i) = 0;
%         FR_area(i) = 0;
%         FR_ave(i) = 0;
%         continue
%     end    
%     a = find(stim_diff==min(stim_diff));
%     FR_temp(i) = max(FR_dynamic(event_start(a):event_end(a)));
%     FR_area(i) = sum(FR_dynamic(event_start(a):event_end(a)))*0.001;
%     FR_ave(i) = mean(FR_dynamic(event_start(a):event_end(a)));
%     event_start_stim(i) = round(T_down(event_start(a)))/dt;
%     event_end_stim(i) = round(T_down(event_end(a))/dt);
% end
% num = zeros(1,n_stim);
% precision = [];
% spike_within_event = {};
% for j = 1:length(event_start_stim)
%     t_event = [];
%     for i = 1:n_recording
%         spike_within_event{i,j} = find(spike_point{i}>=event_start_stim(j) & spike_point{i}<=event_end_stim(j));
%         if ~isempty(spike_within_event{i,j})
%             num(j) = num(j)+1;
%             t_event = [t_event; dt*spike_within_event{i,j}];
%         end
%     end
%     precision(j) = std(t_event);
% end
% reliability = num/10;
% p_spike = zeros(1,size(spike_time_trial,2));
% for i = 1:size(spike_time_trial,2)
% for j = 1:size(spike_time_trial,1)
% if ~isempty(spike_time_trial{j,i})
% p_spike(i) = p_spike(i)+1;
% end
% end
% end
% figure
% plot(T_down, FR_dynamic)
% p_spike = p_spike/10;

p_spike = zeros(1, n_stim);
for i = 1:size(spike_time_trial,2)
    for j = 1:size(spike_time_trial,1)
        if ~isempty(spike_time_trial{j,i})
            p_spike(i) = p_spike(i)+1;
        end
    end
end
p_spike = p_spike/size(spike_time_trial,1);

jitter = zeros(1, n_stim);
for i = 1:size(spike_time_trial,2)
    a = [];
    for j = 1:size(spike_time_trial,1)
        if ~isempty(spike_time_trial{j,i})
            a = [a, min(spike_time_trial{j,i})*1e3];
        end
        if ~isempty(a)
            jitter(i) = std(a);
        end
    end
end
p_burst = zeros(1, n_stim);
for i = 1:size(spike_time_trial,2)
    for j = 1:size(spike_time_trial,1)
        if length(spike_time_trial{j,i})>=2
            p_burst(i) = p_burst(i)+1;
        end
    end
end
p_burst = p_burst/size(spike_time_trial,1);
save(fullfile(path, sprintf('summation_rep_new_%d.mat', idx_SLM)),'jitter','Vm','Vm_seg','Im', 'T','idx_SLM','dt', 'stim_start', 'stim_start_time', 'spike_time', 'FR', 'spike_time_trial', 'mean_FR','p_spike', 'p_burst')
