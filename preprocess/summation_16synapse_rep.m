%%
close all
% clear all

idx_single =[10];
spine = {[1:16]};
idx_sum = 10;
idx_SLM = 1;
branch_no =1;


path = 'E:\data\uncaging\multi cluster cooperation\013123\cell2';
preflix = '230131_001';

clear Vm PC shutter
filename_single = sprintf('%s.single_spine_series_16synapses.%d.wcp',preflix, idx_single(1));
%filename_sum_SLM = sprintf('%s.summation.%d.wcp',preflix, idx_sum_SLM);

out=import_wcp(fullfile(path, filename_single),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
for i = 1
    Vm(:,i) = out.S{3}(:,i);
    PC(:,i) = out.S{7}(:,i);
    stim_start{i} = find(diff(PC(:,i))>90);
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>40)+1]);
    if size(stim_start{i}, 1) == 0
        continue
    end
    num_spines = length(stim_start{i});       
end
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP = [];
if length(idx_single) == 1
figure
for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>50);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>50)+1]);
    for j = 1:length(stim_start{i})
        if shutter(stim_start{i}(j), i)< 4
            continue
        else
            stim_start_new = 500;
            EPSP(:,j) = Vm(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i)-Vm(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp(j),max_amp(j)] = max(abs(EPSP(stim_start_new:stim_start_new+500,j)-EPSP(stim_start_new,j)));
            max_amp(j) = max_amp(j)+stim_start_new;
            [half,I] = sort(abs(EPSP(:,j)-amp(j)/2));
            b = find(I<max_amp(j) & I>stim_start_new);
            a = find(I>max_amp(j));
            halfwidth(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime(j) = (max_amp(j)-stim_start_new)*dt*1000; %in ms
            maxdv(j) = max(diff(movmean(EPSP(:,j), 10)))/dt/1000; %in V/s
            subplot(4,4,j)
            plot(t, EPSP(:,j), 'Color', [0.5,0.5,0.5])
                xlim([-20,150])
%             ylim([-0.5,1.5])
            hold on
%             scatter(t([I(a(1)),I(b(1)),max_amp(j)]), EPSP([I(a(1)),I(b(1)),max_amp(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])            
            xlabel('t (ms)')
            ylabel('Vm (mV)')
%             xlim([min(t), max(t)])
            box off
        end
    end
end
else
    for m = 1:length(idx_single)
        filename_single = sprintf('%s.single_spine_series.%d.wcp',preflix, idx_single(m));
        out=import_wcp(fullfile(path, filename_single),'debug');
        for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    PC(:,i) = out.S{7}(:,i);
    shutter(:,i) = out.S{8}(:,i);
    stim_start{i} = find(diff(PC(:,i))>50);
    if size(stim_start{i}, 1) == 0
        continue
    end
    stim_start{i} = stim_start{i}([1; find(diff(stim_start{i})>50)+1]);
    for j = spine{m}%1:length(stim_start{i})
        if shutter(stim_start{i}(j), i)< 4
            continue
        else
            stim_start_new = 500;
            EPSP(:,j) = Vm(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4001,i)-Vm(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp(j),max_amp(j)] = max(abs(EPSP(stim_start_new:stim_start_new+800,j)-EPSP(stim_start_new,j)));
            max_amp(j) = max_amp(j)+stim_start_new;
            [half,I] = sort(abs(EPSP(:,j)-amp(j)/2));
            b = find(I<max_amp(j) & I>stim_start_new);
            a = find(I>max_amp(j));
            halfwidth(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime(j) = (max_amp(j)-stim_start_new)*dt*1000; %in ms
            maxdv(j) = max(diff(movmean(EPSP(:,j), 10)))/dt/1000; %in V/s
            subplot(2,5,j)
            plot(t, EPSP(:,j), 'Color', [0.5,0.5,0.5])
                        xlim([-20,150])
            ylim([-0.5,1.5])
            hold on
%             scatter(t([I(a(1)),I(b(1)),max_amp(j)]), EPSP([I(a(1)),I(b(1)),max_amp(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
            xlabel('t (ms)')
            ylabel('Vm (mV)')
%             xlim([min(t), max(t)])
            box off
        end
    end
end
        
    end
end


save(fullfile(path, sprintf('single_spine_series_16syn_%d.mat', idx_single(1))),'EPSP', 'Vm', 'amp', 'risetime','halfwidth', 'maxdv','idx_single', 't', 'dt')

%% 16synapse_gradual_sum
path = 'E:\data\uncaging\multi cluster cooperation\011823\cell2';
preflix = '230118_001';
close all
idx_SLM = 37;
clear PC shutter Vm_SLM Im_SLM EPSP_SLM
filename_SLM = sprintf('%s.summation_16synapses.%d.wcp',preflix, idx_SLM);
out=import_wcp(fullfile(path, filename_SLM),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP_SLM = [];
figure
for i = 1:n_recording
    Vm_SLM(:,i) = out.S{3}(:,i);
    Im_SLM(:,i)  = out.S{4}(:,i);
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
            EPSP_SLM(:,j) = Vm_SLM(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i);%-Vm_SLM(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp_SLM(j),max_amp_SLM(j)] = max(abs(EPSP_SLM(stim_start_new:stim_start_new+1200,j)-EPSP_SLM(stim_start_new,j)));
            max_amp_SLM(j) = max_amp_SLM(j)+stim_start_new;
            [half,I] = sort(abs(EPSP_SLM(:,j)-amp_SLM(j)/2));
            b = find(I<max_amp_SLM(j) & I>stim_start_new);
            a = find(I>max_amp_SLM(j));
            halfwidth_SLM(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_SLM(j) = (max_amp_SLM(j)-stim_start_new)*dt*1000; %in ms
            maxdv_SLM(j) = max(diff(movmean(EPSP_SLM(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_SLM(:,j), 10)));
            subplot(3,5,j)
            plot(t, EPSP_SLM(:,j), 'Color', [0.5,0.5,0.5])
            hold on
            scatter(t([I(a(1)),I(b(1)),max_amp_SLM(j)]), EPSP_SLM([I(a(1)),I(b(1)),max_amp_SLM(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
            scatter(t(c), EPSP_SLM(c, j))
            xlabel('t (ms)')
            ylabel('Vm (mV)')
            xlim([min(t), max(t)])
            box off
        end
    end
end
% Vm_SLM(stim_start{1}(1), 1)
save(fullfile(path, sprintf('sum_16syn_%d.mat', idx_SLM)),'EPSP_SLM', 'Vm_SLM', 'amp_SLM', 'risetime_SLM','halfwidth_SLM', 'maxdv_SLM','t','idx_SLM','dt')

%%
close all
% clear all


idx_SLM = 1;
branch_no = 1;


path = 'E:\data\uncaging\multi cluster cooperation\110622\cell4';
preflix = '221107_001';


clear PC shutter
filename_SLM = sprintf('%s.summation_16synapses_rep.%d.wcp',preflix, idx_SLM);
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
n_stim  = 15;
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
spike_time = {};
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

window_length = 0.1;
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
save(fullfile(path, sprintf('summation_16synapses_rep_%d.mat', idx_SLM)),'jitter','Vm','Vm_seg','Im', 'T','idx_SLM','dt', 'stim_start', 'stim_start_time', 'spike_time', 'FR', 'spike_time_trial', 'mean_FR','p_spike', 'p_burst')
