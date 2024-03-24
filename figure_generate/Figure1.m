%% Figure 1Aii & S1A
load(fullfile(pwd, 'stat', 'Figure1Aii.mat'))
figure
scatter(galvo_x_dist, galvo_x_raw,'o', 'k')
hold on
scatter(galvo_z_dist, galvo_z_raw,'o', 'r')
plot(galvo_x_dist_fit, galvo_x_fit,'Color', 'k')
plot(galvo_z_dist_fit, galvo_z_fit,'Color','r')
xlabel('Distance (um)')
ylabel('Intensity [A.U.]')
ylim([0,1.1])
xlim([-5,5])
box on
title('galvo PSF')

figure
scatter(SLM_x_dist, SLM_x_raw,'o', 'k')
hold on
scatter(SLM_z_dist, SLM_z_raw,'o', 'r')
plot(SLM_x_dist_fit, SLM_x_fit,'Color', 'k')
plot(SLM_z_dist_fit, SLM_z_fit,'Color','r')
xlabel('Distance (um)')
ylabel('Intensity [A.U.]')
ylim([0,1.1])
xlim([-5,5])
box on
title('SLM PSF')
%% Figure 1Ci
load(fullfile(pwd, 'stat', 'Figure1Ci.mat'))
figure
subplot(121)
scatter(offset_lateral, amp_lateral, 'k')
hold on
xplot_lateral = 0:0.01:3;
plot(xplot_lateral, f_SLM_lateral(xplot_lateral), 'k')
ylim([0,1.1])
xlim([0,3])
subplot(122)
scatter(offset_axial, amp_axial, 'r')
hold on
xplot_axial = -5:0.01:5;
plot(xplot_axial, f_SLM_axial(xplot_axial), 'r')
ylim([0,1.1])
xlim([-5,5])

%% Figure 1Cii
load(fullfile(pwd, 'stat', 'Figure1Cii.mat'))
figure
scatter(power_g_g_sort,amp_g_g_sort, 'o', 'b')
xplot1 = power_g_g_sort(1):0.01:power_g_g_sort(end);
hold on
plot(xplot1, exp(log(xplot1)*P_g_g(1)+P_g_g(2)), 'b')
scatter(power_com_sort,amp_com_sort, 'o', 'r')
xplot1 = power_com_sort(1):0.01:power_com_sort(end);
hold on
plot(xplot1, exp(log(xplot1)*P_com(1)+P_com(2)), 'r')
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
box on
%% Figure 1Diii
load(fullfile(pwd, 'stat', 'Figure1D.mat'))
offset = offset(6:20);
amps = amps(:,6:20);
curve = fit(offset', mean(amps, 1)', 'poly4');
addpath(genpath('plotting'));
errorbar_with_fitcurve(offset, {amps'}, curve);
xlim([-40, 40])
%% Figure 1E
load(fullfile(pwd, 'stat', 'Figure1E.mat'))
figure, plot(t, EPSP, 'k')
figure, plot(t_Ca, raw_df, 'Color', [0.5,0.5,0.5]), hold on
plot(t_Ca, filt_df, 'Linewidth', 2)
%% Figure S1D
load(fullfile(pwd, 'stat', 'FigureS1D.mat'))
% clear all
% EPSP_SLM = {};
% EPSP_g_g = {};
% dt = [];
% 
% % Iclamp
% date = [100920, 100220, 100220, 100220, 092820, 092820, 092820, 092820, 092620, 092620, 092620, 092620,092620, 092620, 092420, 091320, 091320];
% cell_idx = [2,2,1,1,2,2,1,1,3,3,3,2,2,2,2,4,3];
% branch_idx = [3,1,1,3,1,3,1,2,6,1,2,2,8,9,2,1,1];
% 
% %Vclamp
% % date = [100920, 100920,100920,100220,100220,092820,092820,092820,092420,091920,090720];
% % cell_idx = [3,2,1,2,1,2,2,2,2,4,4];
% % branch_idx = [2,4,1,3,2,4,5,6,2,1,1];
% 
% for i = 1:length(date)
% path = sprintf('Z:\\data\\shulan\\uncaging\\physiological characterization\\%06d\\cell%d', date(i), cell_idx(i));
% filename_SLM = sprintf('branch%d_SLM.mat', branch_idx(i));
% filename_g_g = dir(fullfile(path,sprintf('branch%d_g-g*.mat', branch_idx(i))));
% load(fullfile(path,filename_g_g.name));
% EPSP_SLM{i} = EPSP_2spine;
% load(fullfile(path, filename_SLM));
% EPSP_g_g{i} = EPSP_2spine;
% dt = [dt, t(2)-t(1)];
% end
% 
% 
% for i = 1:length(date)
%     stim_start = 500;%find(EPSP_g_g{i}==0);
%     amp_g_g(i) = max(abs(EPSP_g_g{i}(stim_start:stim_start+500))-EPSP_g_g{i}(stim_start));
%     rise_half = find(abs(EPSP_g_g{i}-amp_g_g(i)/2)<=5e-2);
%     b = find(rise_half<stim_start+200);
%     fall_half = find(abs(EPSP_g_g{i}-amp_g_g(i)/2)<=1e-2);
%     max_amp = find(EPSP_g_g{i}==max(abs(EPSP_g_g{i}(stim_start:stim_start+500))));
%     a = find(fall_half>stim_start+400);
%     halfwidth_g_g(i) = (fall_half(a(1))-rise_half(b(end)))*dt(i);
%     risetime_g_g(i) = (max_amp(1)-stim_start)*dt(i);
%     
%     amp_SLM(i) = max(abs(EPSP_SLM{i}(stim_start:stim_start+500))-EPSP_SLM{i}(stim_start));
%     rise_half = find(abs(EPSP_SLM{i}-amp_SLM(i)/2)<=5e-2);
%     b = find(rise_half<stim_start+200);
%     fall_half = find(abs(EPSP_SLM{i}-amp_SLM(i)/2)<=1e-2);
%     max_amp = find(EPSP_SLM{i}==max(abs(EPSP_SLM{i}(stim_start:stim_start+500))));
%     a = find(fall_half>stim_start+400);
%     halfwidth_SLM(i) = (fall_half(a(1))-rise_half(b(end)))*dt(i);
%     risetime_SLM(i) = (max_amp(1)-stim_start)*dt(i);
% end

addpath(genpath('plotting'));
boxplot_compact({amp_SLM, amp_g_g}, ['r', 'k'])
boxplot_compact({halfwidth_SLM, halfwidth_g_g}, ['r', 'k'])


% save(fullfile(pwd, 'stat','FigureS1D.mat'), 'EPSP_g_g', 'EPSP_SLM', 'dt', 'amp_g_g', 'amp_SLM', 'risetime_g_g', 'risetime_SLM', 'halfwidth_g_g', 'halfwidth_SLM')
%% Figure S1E
load(fullfile(pwd, 'stat', 'FigureS1E.mat'))
% clear all
% EPSP_SLM = {};
% EPSP_g_g = {};
% dt = [];
% 
% % % Iclamp
% % date = [100920, 100220, 100220, 100220, 092820, 092820, 092820, 092820, 092620, 092620, 092620, 092620,092620, 092620, 092420, 091320, 091320];
% % cell_idx = [2,2,1,1,2,2,1,1,3,3,3,2,2,2,2,4,3];
% % branch_idx = [3,1,1,3,1,3,1,2,6,1,2,2,8,9,2,1,1];
% 
% % Vclamp
% date = [100920, 100920,100920,100220,100220,092820,092820,092820,092420,091920,090720];
% cell_idx = [3,2,1,2,1,2,2,2,2,4,4];
% branch_idx = [2,4,1,3,2,4,5,6,2,1,1];
% 
% for i = 1:length(date)
% path = sprintf('Z:\\data\\shulan\\uncaging\\physiological characterization\\%06d\\cell%d', date(i), cell_idx(i));
% filename_SLM = sprintf('branch%d_SLM.mat', branch_idx(i));
% filename_g_g = dir(fullfile(path,sprintf('branch%d_g-g*.mat', branch_idx(i))));
% load(fullfile(path,filename_g_g.name));
% EPSP_SLM{i} = EPSP_2spine;
% load(fullfile(path, filename_SLM));
% EPSP_g_g{i} = EPSP_2spine;
% dt = [dt, t(2)-t(1)];
% end
% 
% 
% for i = 1:length(date)
%     stim_start = 500;%find(EPSP_g_g{i}==0);
%     amp_g_g(i) = max(abs(EPSP_g_g{i}(stim_start:stim_start+500))-EPSP_g_g{i}(stim_start));
%     rise_half = find(abs(EPSP_g_g{i}-amp_g_g(i)/2)<=5e-2);
%     b = find(rise_half<stim_start+200);
%     fall_half = find(abs(EPSP_g_g{i}-amp_g_g(i)/2)<=1e-2);
%     max_amp = find(EPSP_g_g{i}==max(abs(EPSP_g_g{i}(stim_start:stim_start+500))));
%     a = find(fall_half>stim_start+400);
%     halfwidth_g_g(i) = (fall_half(a(1))-rise_half(b(end)))*dt(i);
%     risetime_g_g(i) = (max_amp(1)-stim_start)*dt(i);
%     
%     amp_SLM(i) = max(abs(EPSP_SLM{i}(stim_start:stim_start+500))-EPSP_SLM{i}(stim_start));
%     rise_half = find(abs(EPSP_SLM{i}-amp_SLM(i)/2)<=5e-2);
%     b = find(rise_half<stim_start+200);
%     fall_half = find(abs(EPSP_SLM{i}-amp_SLM(i)/2)<=1e-2);
%     max_amp = find(EPSP_SLM{i}==max(abs(EPSP_SLM{i}(stim_start:stim_start+500))));
%     a = find(fall_half>stim_start+400);
%     halfwidth_SLM(i) = (fall_half(a(1))-rise_half(b(end)))*dt(i);
%     risetime_SLM(i) = (max_amp(1)-stim_start)*dt(i);
% end

boxplot_compact({amp_SLM, amp_g_g}, ['r', 'k'])
boxplot_compact({halfwidth_SLM, halfwidth_g_g}, ['r', 'k'])

% save(fullfile(pwd, 'stat','FigureS1D.mat'), 'EPSP_g_g', 'EPSP_SLM', 'dt', 'amp_g_g', 'amp_SLM', 'risetime_g_g', 'risetime_SLM', 'halfwidth_g_g', 'halfwidth_SLM')
