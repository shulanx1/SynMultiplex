close all
clear all

idx_single =[2];
spine = {[1:10]};
idx_sum = 11;
idx_SLM = 1;
branch_no =1;


path = 'E:\data\dendritic patch\082723\cell3';
preflix = '230827_001';
%%
filename_single = sprintf('%s.single_spine_series.%d.wcp',preflix, idx_single(1));
%filename_sum_SLM = sprintf('%s.summation.%d.wcp',preflix, idx_sum_SLM);

clear Vm Im PC shutter
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
            subplot(2,4,j)
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


save(fullfile(path, sprintf('single_spine_series_%d.mat', idx_single(1))),'EPSP', 'Vm', 'amp', 'risetime','halfwidth', 'maxdv','idx_single', 't', 'dt')
%%
clear PC shutter
filename_sum = sprintf('%s.summation.%d.wcp',preflix, idx_sum);
out=import_wcp(fullfile(path, filename_sum),'debug');
n_channel = out.channel_no;
n_recording = length(out.rec_index);

stim_start = cell(1,n_recording);
T = out.T;
dt = T(2) - T(1);
t = (-499*dt:dt:4501*dt)*1e3;

EPSP_sum = [];
figure
for i = 1:n_recording
    Vm_sum(:,i) = out.S{3}(:,i);
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
            EPSP_sum(:,j) = Vm_sum(stim_start{i}(j)-stim_start_new+1:stim_start{i}(j)+4501,i)-Vm_sum(stim_start{i}(j), i);
            %EPSP = [EPSP, Vm(stim_start{i}(j)-500:stim_start{i}(j)+3500,i)];
            [amp_sum(j),max_amp_sum(j)] = max(abs(EPSP_sum(stim_start_new:stim_start_new+1200,j)-EPSP_sum(stim_start_new,j)));
            max_amp_sum(j) = max_amp_sum(j)+stim_start_new;
            [half,I] = sort(abs(EPSP_sum(:,j)-amp_sum(j)/2));
            b = find(I<max_amp_sum(j) & I>stim_start_new);
            a = find(I>max_amp_sum(j));
            halfwidth_sum(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_sum(j) = (max_amp_sum(j)-stim_start_new)*dt*1000; %in ms
            maxdv_sum(j) = max(diff(movmean(EPSP_sum(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_sum(:,j), 10)));
            subplot(2,4,j)
            plot(t, EPSP_sum(:,j), 'Color', [0.5,0.5,0.5])
            hold on
            scatter(t([I(a(1)),I(b(1)),max_amp_sum(j)]), EPSP_sum([I(a(1)),I(b(1)),max_amp_sum(j)],j), 'O', 'MarkerEdgeColor', [0, 0.4470, 0.7410])
            scatter(t(c), EPSP_sum(c, j))
            xlabel('t (ms)')
            ylabel('Vm (mV)')
            xlim([min(t), max(t)])
            box off
        end
    end
end
% maxdv_sum(end)
save(fullfile(path, sprintf('sum_galvo_%d.mat', idx_sum)),'EPSP_sum', 'Vm_sum', 'amp_sum', 'risetime_sum','halfwidth_sum', 'maxdv_sum', 'idx_sum', 't','dt')
%%
idx_SLM = 1;
clear PC shutter Vm_SLM Im_SLM amp_SLM EPSP_SLM halfwidth_SLM maxdv_SLM risetime_SLM
filename_SLM = sprintf('%s.summation.%d.wcp',preflix, idx_SLM);
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
            [half,I] = sort(abs((EPSP_SLM(:,j)-Vm_SLM(stim_start{i}(j), i))-amp_SLM(j)/2));
            b = find(I<max_amp_SLM(j) & I>stim_start_new);
            a = find(I>max_amp_SLM(j));
            halfwidth_SLM(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_SLM(j) = (max_amp_SLM(j)-stim_start_new)*dt*1000; %in ms
            maxdv_SLM(j) = max(diff(movmean(EPSP_SLM(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_SLM(:,j), 10)));
            subplot(2,5,j)
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
save(fullfile(path, sprintf('sum_SLM_%d.mat', idx_SLM)),'EPSP_SLM', 'Vm_SLM', 'amp_SLM', 'risetime_SLM','halfwidth_SLM', 'maxdv_SLM','t','idx_SLM','dt')
%%
idx_SLM = 5;
clear PC shutter Vm_SLM Im_SLM amp_SLM EPSP_SLM halfwidth_SLM maxdv_SLM risetime_SLM
filename_SLM = sprintf('%s.summation_ca.%d.wcp',preflix, idx_SLM);
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
            [half,I] = sort(abs((EPSP_SLM(:,j)-Vm_SLM(stim_start{i}(j), i))-amp_SLM(j)/2));
            b = find(I<max_amp_SLM(j) & I>stim_start_new);
            a = find(I>max_amp_SLM(j));
            halfwidth_SLM(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
            risetime_SLM(j) = (max_amp_SLM(j)-stim_start_new)*dt*1000; %in ms
            maxdv_SLM(j) = max(diff(movmean(EPSP_SLM(:,j), 10)))/dt/1000; %in V/s
            [~,c] = max(diff(movmean(EPSP_SLM(:,j), 10)));
            subplot(2,5,j)
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
save(fullfile(path, sprintf('sum_%d.mat', idx_SLM)),'EPSP_SLM', 'Vm_SLM', 'amp_SLM', 'risetime_SLM','halfwidth_SLM', 'maxdv_SLM','t','idx_SLM','dt')
%%
figure
subplot(132)
len = num_spines;
red = [0.9,0.9,0.9];
pink = [0, 0.4470, 0.7410];
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
for i = 1:num_spines
plot(t, EPSP_sum(:,i), 'LineWidth', 1, 'Color', colors_p(i,:))
hold on
end
y1 = ylim();
subplot(131)
len = num_spines;
red = [0.9,0.9,0.9];
pink = [0,0,0];
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
for i = 1:num_spines
plot(t, sum(EPSP(:,1:i),2), 'LineWidth', 1, 'Color', colors_p(i,:))
hold on
end
ylim(y1)
% y1 = ylim();
subplot(133)
len = num_spines;
red = [0.9,0.9,0.9];
pink = [0.8500, 0.3250, 0.0980];
colors_p = [linspace(red(1),pink(1),len)', linspace(red(2),pink(2),len)', linspace(red(3),pink(3),len)'];
plot(t, EPSP(:,1), 'LineWidth', 1, 'Color', colors_p(i,:))
for i = 1:num_spines-1
plot(t, EPSP_SLM(:,i), 'LineWidth', 1, 'Color', colors_p(i+1,:))
hold on
end
ylim(y1)
saveas(gcf, fullfile(path, sprintf('summation_branch%d.eps', branch_no)), 'epsc')
out_file = fullfile(path, sprintf('summation_branch%d.mat', branch_no));
save(out_file, 'EPSP', 'EPSP_sum', 'EPSP_SLM','idx_single', 'idx_sum', 'idx_SLM','t','dt','amp', 'amp_sum', 'amp_SLM','branch_no')%, 'amp', 'halfwidth', 'risetime');
% save(out_file, 'EPSP','EPSP_SLM','idx_single', 'idx_SLM','t','dt','amp', 'amp_SLM','branch_no')%, 'amp', 'halfwidth', 'risetime');
% save(out_file, 'EPSP','EPSP_sum','idx_single', 'idx_sum','t','dt','amp', 'amp_sum','branch_no')%, 'amp', 'halfwidth', 'risetime');
%% pool
% stim_start_new = 500;
% for i = 1:10
%     for j = 1:size(EPSP_single_pool{i}, 2)
%         amp_single_pool{i}(j) = max(abs(EPSP_single_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_single_pool{i}(stim_start_new,j)));
%     end
%         for j = 1:size(EPSP_galvo_pool{i}, 2)
%         amp_galvo_pool{i}(j) = max(abs(EPSP_galvo_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_galvo_pool{i}(stim_start_new,j)));
%         end
%         for j = 1:size(EPSP_SLM_pool{i}, 2)
%         amp_SLM_pool{i}(j) = max(abs(EPSP_SLM_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_SLM_pool{i}(stim_start_new,j)));
%         end
% end
% for i = 16:54
%     amp_SLM_pool{i} = [amp_single_pool{i}(1), amp_SLM_pool{i}];
% end

for i = 1:54
    for j = 1:size(EPSP_single_pool{i}, 2)
        [amp_single_pool{i}(j), max_amp] = max(abs(EPSP_single_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_single_pool{i}(stim_start_new,j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs(EPSP_single_pool{i}(:,j)-amp_single_pool{i}(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        halfwidth_single_pool{i}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
    end
    for j = 1:size(EPSP_galvo_pool{i}, 2)
        [amp_galvo_pool{i}(j), max_amp] = max(abs(EPSP_galvo_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_galvo_pool{i}(stim_start_new,j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs(EPSP_galvo_pool{i}(:,j)-amp_galvo_pool{i}(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        halfwidth_galvo_pool{i}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
    end
    for j = 1:size(EPSP_SLM_pool{i}, 2)
        [amp_SLM_pool{i}(j), max_amp] = max(abs(EPSP_SLM_pool{i}(stim_start_new:stim_start_new+500,j)-EPSP_SLM_pool{i}(stim_start_new,j)));
        max_amp = max_amp+stim_start_new;
        [half,I] = sort(abs(EPSP_SLM_pool{i}(:,j)-amp_SLM_pool{i}(j)/2));
        b = find(I<max_amp & I>stim_start_new);
        a = find(I>max_amp);
        halfwidth_SLM_pool{i}(j) = (I(a(1))-I(b(1)))*dt*1000; %in ms
    end
end
%%
ft = fittype('a+(b-a)/(1+exp(-(x-c)/d))', 'indep', 'x');
curve_galvo_dist = {};
curve_galvo_prox = {};
figure
for i = 1:54
    if isempty(amp_galvo_pool{i})
        continue;
    end
    if dist(i)>=50
        curve_galvo_dist{length(curve_galvo_dist)+1} = fit([1:spine_num(i)]', amp_galvo_pool{i}(1:spine_num(i))',ft,'start', [0,5,4,1]);
        plot(1:spine_num(i), amp_galvo_pool{i}(1:spine_num(i)), 'Color', [0.5,0,0])
        hold on
    else 
        curve_galvo_prox{length(curve_galvo_prox)+1} = fit([1:spine_num(i)]', amp_galvo_pool{i}(1:spine_num(i))',ft,'start', [0,5,4,1]);
        plot(1:spine_num(i), amp_galvo_pool{i}(1:spine_num(i)), 'Color', [0,0,0.5])
        hold on
    end
end
curve_galvo_dist{11} = fit([1:spine_num(15)]', amp_galvo_pool{15}(1:spine_num(15))',ft,'start', [0,5,4,1],'Exclude', 4);
%%
ft = fittype('a+(b-a)/(1+exp(-(x-c)/d))', 'indep', 'x');
curve_SLM_dist = {};
curve_SLM_prox = {};
figure
for i = 1:54
    if length(amp_SLM_pool{i})<=2
        continue;
    end
    if dist(i)>=90
        curve_SLM_dist{length(curve_SLM_dist)+1} = fit([1:spine_num(i)]', amp_SLM_pool{i}(1:spine_num(i))',ft,'start', [30,80,4,1]);
        plot(1:spine_num(i), amp_SLM_pool{i}(1:spine_num(i)), 'Color', [0.5,0,0])
        hold on
    else 
        curve_SLM_prox{length(curve_SLM_prox)+1} = fit([1:spine_num(i)]', amp_SLM_pool{i}(1:spine_num(i))',ft,'start', [30,80,4,1]);
        plot(1:spine_num(i), amp_SLM_pool{i}(1:spine_num(i)), 'Color', [0,0,0.5])
        hold on
    end
end
%%
idx_prox = find(dist<=60);
idx_dist = find((dist>60));
sum_prox = [];
for i = 1:length(idx_prox)
if length(amp_SLM_pool{idx_prox(i)})<7
continue
end
sum_prox_1= [];
for j = 1:7
sum_prox_1(j) = amp_SLM_pool{idx_prox(i)}(j)/mean(amp_single_pool{idx_prox(i)}(1:j));
end
sum_prox = [sum_prox;sum_prox_1];
end


sum_dist = [];
for i = 1:length(idx_dist)
if length(amp_SLM_pool{idx_dist(i)})<7
continue
end
sum_dist_1 = [];
for j = 1:7
sum_dist_1(j) = amp_SLM_pool{idx_dist(i)}(j)/mean(amp_single_pool{idx_dist(i)}(1:j));
end
sum_dist = [sum_dist;sum_dist_1];
end
figure
errorbar([1:7], mean(sum_prox), std(sum_prox)./sqrt(size(sum_prox,1)), 'o','Color', 'k',  'MarkerSize',10,  'MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1])
hold on
errorbar([1:7], mean(sum_dist), std(sum_dist)./sqrt(size(sum_dist,1)), 'o','Color', 'r',  'MarkerSize',10,  'MarkerEdgeColor', 'r','MarkerFaceColor', [1,1,1])

ft = fittype('1+b/(1+exp(-(x-a)/k))', 'independent', 'x');
curve_prox = fit([1:7]',mean(sum_prox)', ft);
curve_dist = fit([1:7]',mean(sum_dist)', ft,'Start', [5,7,2])
x_fit = 1:0.01:7;
hold on
plot(x_fit, curve_prox(x_fit), 'k')
plot(x_fit, curve_dist(x_fit), 'r')