%% Figure 2A
load(fullfile(pwd, 'datasets', 'singleclust.mat'))
addpath(genpath('plotting'))
dt = 5e-5;
i_plot = [137:139];
for i = 1:length(i_plot)
    t = -500*dt:dt:dt*(size(data(i_plot(i)).EPSP, 1)-1-500);
    plot_traces_with_gradient(t, data(i_plot(i)).EPSP);
end
for i = 1:length(i_plot)
    t = -500*dt:dt:dt*(size(data(i_plot(i)).EPSP, 1)-2-500);
    plot_traces_with_gradient(t, diff(data(i_plot(i)).EPSP));
end
%% Figure 2B&C
load(fullfile(pwd, 'stat', 'Figure2B&C.mat'))
idx_strong = find(M.MaxDv_dt>2);
idx_weak = setdiff(1:size(M,1), idx_strong);
figure
scatter(M.Amp(idx_strong), M.MaxDv_dt(idx_strong), 'r')
hold on
scatter(M.Amp(idx_weak), M.MaxDv_dt(idx_weak), 'k')
box on

figure
scatter(M.distToSoma(idx_strong), M.MaxDv_dt(idx_strong), 'r')
hold on
scatter(M.distToSoma(idx_weak), M.MaxDv_dt(idx_weak), 'k')
box on
%% Figure 2D&E
load(fullfile(pwd, 'datasets', 'clusterRMPchange.mat'))
dt = 5e-5;
i_plot = [24, 25, 47];
for i = 1:length(i_plot)
    figure
    subplot(121)
    t = -500*dt:dt:dt*(size(data(i_plot(i)).EPSP, 1)-1-500);
    plot(t, data(i_plot(i)).EPSP, 'k');
    xlim([-0.02,0.2])
    subplot(122)
    t = -500*dt:dt:dt*(size(data(i_plot(i)).stim_response, 1)-1-500);
    plot(t, data(i_plot(i)).stim_response, 'k');
    xlim([-0.02,0.2])
end

boxplot_compact({datapool.jitter_strong/1e3, datapool.jitter_weak/1e3}, ['r', 'k']) % ms
boxplot_compact({datapool.onset_strong/1e3, datapool.onset_weak/1e3}, ['r', 'k']) % ms
boxplot_compact({datapool.area_under_curve_burst, datapool.area_under_curve_single}, [[0,147,69]; [27,117,187]]/256)
%% Figure 2G&H
load(fullfile(pwd, 'stat', 'Figure2M.mat'))
i_plot = [35,36];
dt = 5e-5;
figure
subplot(2,2,[1,3])
t = -500*dt:dt:dt*(size(data(i_plot(1)).EPSP_singleclust, 1)-1-500);
plot(t, data(i_plot(1)).EPSP_singleclust(:,1), 'k')
hold on
plot(t, data(i_plot(1)).EPSP_singleclust(:,2))
plot(t, data(i_plot(2)).EPSP_singleclust(:,2))
xlim([-0.02,0.2])
subplot(222)
t = -500*dt:dt:dt*(size(data(i_plot(1)).EPSP_highcond, 1)-1-500);
plot(t, data(i_plot(1)).EPSP_highcond, 'k')
xlim([-0.01,0.05])
subplot(224)
plot(t, data(i_plot(2)).EPSP_highcond, 'k')
xlim([-0.01,0.05])
%% Figure 2K&L
load(fullfile(pwd, 'stat', 'Figure2M.mat'))
i_plot = [2,3];
dt = 5e-5;
figure
subplot(2,2,[1,3])
t = -500*dt:dt:dt*(size(data(i_plot(1)).EPSP_singleclust, 1)-1-500);
plot(t, data(i_plot(1)).EPSP_singleclust(:,1), 'k')
hold on
plot(t, data(i_plot(1)).EPSP_singleclust(:,2))
plot(t, data(i_plot(2)).EPSP_singleclust(:,2))
xlim([-0.02,0.2])
subplot(222)
t = -500*dt:dt:dt*(size(data(i_plot(1)).EPSP_highcond, 1)-1-500);
plot(t, data(i_plot(1)).EPSP_highcond(:,1), 'k')
xlim([-0.01,0.1])
subplot(224)
plot(t, data(i_plot(2)).EPSP_highcond(:,1), 'k')
xlim([-0.01,0.1])
%% Figure 2M
load(fullfile(pwd, 'stat', 'Figure2M.mat'))
addpath(genpath('plotting'))
boxplot_compact({datapool.jitter_strong/1e3, datapool.jitter_weak/1e3}, ['r', 'k']) % ms
boxplot_compact({datapool.onset_strong/1e3, datapool.onset_weak/1e3}, ['r', 'k']) % ms
boxplot_compact({datapool.area_burst, datapool.area_single}, [[0,147,69]; [27,117,187]]/256)

%% Figure S2A&B
load(fullfile(pwd, 'stat', 'FigureS2AB.mat'))
addpath(genpath('plotting'))
errorbar_with_fitcurve([1:8],{amp_NMDA'}),title('FigureS2A')
errorbar_with_fitcurve([1:8],{amp_Na'}),title('FigureS2B')
%% Figure S2D&E
load(fullfile(pwd, 'datasets', 'singleclust_AP5.mat'))
% load(fullfile(pwd, 'stat', 'FigureS2E.mat'))
addpath(genpath('plotting'))
boxplot_pairwise([datapool.area_con_clust,datapool.area_ap5_clust],[[0,0,0]; [0,0,1]])
boxplot_pairwise([datapool.width_con_clust,datapool.width_ap5_clust],[[0,0,0]; [0,0,1]])
errorbar_with_fitcurve([2:8],{datapool.amp_con_norm',datapool.amp_ap5_norm'}, {datapool.curve_con,datapool.curve_ap5}, [[0,0,0]; [0,0,1]]), box on
%% Figure S2G&H
load(fullfile(pwd, 'datasets', 'singleclust_ttx.mat'))
% load(fullfile(pwd, 'stat', 'FigureS2H.mat'))
addpath(genpath('plotting'))
boxplot_pairwise({[datapool.dv_con_clust',datapool.dv_ttx_clust'],[datapool.dv_con_clust_weak',datapool.dv_ttx_clust_weak']},[[1,0,0]; [0,0,0];[1,0.5,0.5];[0,0,0]])
boxplot_pairwise({[datapool.amp_con_clust',datapool.amp_ttx_clust'],[datapool.amp_con_clust_weak',datapool.amp_ttx_clust_weak']},[[1,0,0]; [0,0,0];[1,0.5,0.5];[0,0,0]])
errorbar_with_fitcurve([2:8],{datapool.amp_con_norm',datapool.amp_ttx_norm', datapool.amp_ap5_norm'}, {datapool.curve_con,datapool.curve_ttx, datapool.curve_ap5}, [[1,0,0];[0,0,0]; [0,0,1]]), box on
%% Figure S3C&G
load(fullfile(pwd, 'stat', 'FigureS3.mat'))
errorbar_with_fitcurve([1:8],{datapool.amp_nomod_norm',datapool.amp_nomod_norm_1'}, {datapool.curve_nomod,datapool.curve_nomod_1}, [[0.5,0.5,0.5];[0,0,0]]), box on, title('FigureS3C')
errorbar_with_fitcurve([1:8],{datapool.amp_nomod_norm',datapool.amp_mod_norm'}, {datapool.curve_nomod,datapool.curve_mod}, [[0,0,0];[1,0,0]]), box on, title('FigureS3G')
%% Figure S4A
load(fullfile(pwd, 'stat', 'FigureS4A.mat'))
boxplot_compact({area_precise, area_jitter}, ['r', 'k']) % ms
%% Figure S4C
% simulation done with python
load(fullfile(pwd, 'stat', 'FigureS4C.mat'))
figure,subplot(221),plot(t, reshape(vm(1,:,:), 10,[])), xlim([0,200]), subplot(222),plot(t, reshape(vm(2,:,:), 10,[])), xlim([0,200]);
subplot(223),plot(t, [zeros(1, 5000),noise_i(1,:),zeros(1, (length(t)-size(noise_i, 2)-5000))], 'r'), xlim([0,200]), subplot(224), plot(t, [zeros(1, 5000),noise_i(2,:),zeros(1, (length(t)-size(noise_i, 2)-5000))], 'k'), xlim([0,200])

idx_plot = 1; % as an example, plot the phase plane plot of the first trial
mh_Na = reshape(mh(1, idx_plot, 8, :), 1, []);
mh_NMDA = reshape(mh(2, idx_plot, 8, :), 1, []);
n_Na = reshape(n(1, idx_plot, 5, :), 1, []);
n_NMDA = reshape(n(2, idx_plot, 5, :), 1, []);

gamma_K = 20e-9; gamma_Na = 20e-9; D_Na = 60e8; D_K = 18e8; NNa = 12000; NK = 3600; % parameters used in the simulation
gna_Na = gamma_Na/(NNa/D_Na)*double(mh_Na);
gk_Na = gamma_K/(NK/D_K)*double(n_Na);

gna_NMDA = gamma_Na/(NNa/D_Na)*double(mh_NMDA);
gk_NMDA = gamma_K/(NK/D_K)*double(n_NMDA);

figure,subplot(121), plot3(gna_Na,gk_Na,reshape(vm(1,idx_plot,:), 1,[]), 'r'),subplot(122), plot3(gna_NMDA,gk_NMDA,reshape(vm(2,idx_plot,:), 1,[]), 'k')
%% Figure S5A-E
load(fullfile(pwd, 'datasets', 'bkgandclust_precision_and_burst.mat'))
i_plot = 43; % figure S5B, for figure S5D, change i_plot tp 29;
dt = 5e-5;
t = -500*dt:dt:dt*(size(data(i_plot).EPSP_singleclust, 1)-1-500);
figure,plot(t, data(i_plot).EPSP_singleclust), xlim([-0.02,0.2]);
t = -500*dt:dt:dt*(size(data(i_plot).EPSP_highcond, 1)-1-500);
figure,plot(t, data(i_plot).EPSP_highcond, 'k'), xlim([-0.02,0.2]);

maxdv = [data.maxdv_single];
maxdv_clust = maxdv(2:2:end);
boxplot_compact({maxdv_clust(datapool.idx_burst), maxdv_clust(datapool.idx_single)}, [[0,147,69]; [27,117,187]]/256)
hold on
for i = 1:length(datapool.idx_burst)
    if maxdv_clust(datapool.idx_burst(i))>2
        scatter(1.2, maxdv_clust(datapool.idx_burst(i)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None')
    else
        scatter(1.2, maxdv_clust(datapool.idx_burst(i)), 'MarkerFaceColor', [0.5,0.5,0.5], 'MarkerEdgeColor', 'None')
    end
end
for i = 1:length(datapool.idx_single)
    if maxdv_clust(datapool.idx_single(i))>2
        scatter(2.2, maxdv_clust(datapool.idx_single(i)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'None')
    else
        scatter(2.2, maxdv_clust(datapool.idx_single(i)), 'MarkerFaceColor', [0.5,0.5,0.5], 'MarkerEdgeColor', 'None')
    end
end

%% Figure S5G
load(fullfile(pwd, 'stat', 'FigureS5G.mat'))
boxplot_with_datapoint(amp_temp(:,2:3))
boxplot_with_datapoint(hw_temp(:,2:3))
%% Figure S6
load(fullfile(pwd,'datasets','highcond_multibranch_allconditions.mat'));
load(fullfile(pwd,'datasets', 'celltypes.mat'));
label = M.idx;
label1 = zeros(length(label), 1);
label1(label~=3)=1;

a = {};
for i = 1:size(M,1)
a(i) = M.date(i);
end

idx = [];
for i = 1:length(data)
    dup = strcmp(data(i).date, a);
    if any(dup)
        idx_dup = find(dup);
        for j = 1:length(idx_dup)
            if strcmp(data(i).cell, num2str(M.cell(idx_dup(j))))
                data(i).celltype = label1(idx_dup(j));
                if label(idx_dup(j))==2
                    data(i).IB = 1;
                else
                    data(i).IB = 0;
                end
            end
        end
    end
    if isempty(data(i).celltype)
        idx = [idx, i];       
    end
end

celltype = [data.celltype];
IB = [data.IB];
idx_l5a = find(celltype==0);
idx_RS = find((celltype==1)&(IB==0));
idx_IB = find((celltype==1)&(IB==1));
date = {};
for i = 1:length(data)
date{i} = [data(i).date, data(i).cell];
end

p_burst = [data.p_burst];
for i = 1:length(data)
if data(i).p_spike(1)==0
p_burst(i) = data(i).p_burst(2);
else
p_burst(i) = data(i).p_burst(1);
end
end

IB_burst = unique(date(idx_IB(find(p_burst(idx_IB)>0))));
IB_cell = unique(date(idx_IB));
idx_IB_strong = [];
for i = 1:length(idx_IB)
if ~isempty(find(datapool.idx_strong==idx_IB(i)))
idx_IB_strong = [idx_IB_strong, idx_IB(i)];
end
end
IB_strong = unique(date(idx_IB_strong));

RS_burst = unique(date(idx_RS(find(p_burst(idx_RS)>0))));
RS_cell = unique(date(idx_RS));
idx_RS_strong = [];
for i = 1:length(idx_RS)
if ~isempty(find(datapool.idx_strong==idx_RS(i)))
idx_RS_strong = [idx_RS_strong, idx_RS(i)];
end
end
RS_strong = unique(date(idx_RS_strong));

l5a_burst = unique(date(idx_l5a(find(p_burst(idx_l5a)>0))));
l5a_cell = unique(date(idx_l5a));
idx_l5a_strong = [];
for i = 1:length(idx_l5a)
if ~isempty(find(datapool.idx_strong==idx_l5a(i)))
idx_l5a_strong = [idx_l5a_strong, idx_l5a(i)];
end
end
l5a_strong = unique(date(idx_l5a_strong));
%% find cells that generated burst/spikelet under bkg + cluster
load(fullfile(pwd,'datasets','bkgandclust_precision_and_burst.mat'));

a = {};
for i = 1:size(M,1)
a(i) = M.date(i);
end

idx = [];
for i = 1:length(data)
    dup = strcmp(data(i).date, a);
    if any(dup)
        idx_dup = find(dup);
        for j = 1:length(idx_dup)
            if strcmp(data(i).cell, num2str(M.cell(idx_dup(j))))
                data(i).celltype = label1(idx_dup(j));
                if label(idx_dup(j))==2
                    data(i).IB = 1;
                else
                    data(i).IB = 0;
                end
            end
        end
    end
    if isempty(data(i).celltype)
        idx = [idx, i];       
    end
end

celltype = [data.celltype];
IB = [data.IB];
idx_l5a = find(celltype==0);
idx_RS = find((celltype==1)&(IB==0));
idx_IB = find((celltype==1)&(IB==1));
date = {};
for i = 1:length(data)
date{i} = [data(i).date, data(i).cell];
end


p_burst = [data.p_burst];
IB_burst1 = unique(date(idx_IB(find(p_burst(idx_IB)>0))));
IB_cell1 = unique(date(idx_IB));


RS_burst1 = unique(date(idx_RS(find(p_burst(idx_RS)>0))));
RS_cell1 = unique(date(idx_RS));


l5a_burst1 = unique(date(idx_l5a(find(p_burst(idx_l5a)>0))));
l5a_cell1 = unique(date(idx_l5a));


IB_cell1 = unique(date(idx_IB));
idx_IB_strong = [];
for i = 1:length(idx_IB)
if ~isempty(find(datapool.idx_strong==idx_IB(i)))
idx_IB_strong = [idx_IB_strong, idx_IB(i)];
end
end
IB_strong1 = unique(date(idx_IB_strong));


RS_cell1 = unique(date(idx_RS));
idx_RS_strong = [];
for i = 1:length(idx_RS)
if ~isempty(find(datapool.idx_strong==idx_RS(i)))
idx_RS_strong = [idx_RS_strong, idx_RS(i)];
end
end
RS_strong1 = unique(date(idx_RS_strong));


l5a_cell1 = unique(date(idx_l5a));
idx_l5a_strong = [];
for i = 1:length(idx_l5a)
if ~isempty(find(datapool.idx_strong==idx_l5a(i)))
idx_l5a_strong = [idx_l5a_strong, idx_l5a(i)];
end
end
l5a_strong1 = unique(date(idx_l5a_strong));
%% find cells that generated spikelet under single cluster stimulation
load(fullfile(pwd,'datasets','singleclust.mat'));

a = {};
for i = 1:size(M,1)
a(i) = M.date(i);
end

idx = [];
for i = 1:length(data)
    dup = strcmp(data(i).date, a);
    if any(dup)
        idx_dup = find(dup);
        for j = 1:length(idx_dup)
            if strcmp(data(i).cell, num2str(M.cell(idx_dup(j))))
                data(i).celltype = label1(idx_dup(j));
                if label(idx_dup(j))==2
                    data(i).IB = 1;
                else
                    data(i).IB = 0;
                end
            end
        end
    end
    if isempty(data(i).celltype)
        idx = [idx, i];       
    end
end

celltype = [data.celltype];
IB = [data.IB];
idx_l5a = find(celltype==0);
idx_RS = find((celltype==1)&(IB==0));
idx_IB = find((celltype==1)&(IB==1));
date = {};
for i = 1:length(data)
date{i} = [data(i).date, data(i).cell];
end


IB_cell2 = unique(date(idx_IB));
idx_IB_strong = [];
for i = 1:length(idx_IB)
if ~isempty(find(datapool.strong_idx==idx_IB(i)))
idx_IB_strong = [idx_IB_strong, idx_IB(i)];
end
end
IB_strong2 = unique(date(idx_IB_strong));


RS_cell2 = unique(date(idx_RS));
idx_RS_strong = [];
for i = 1:length(idx_RS)
if ~isempty(find(datapool.strong_idx==idx_RS(i)))
idx_RS_strong = [idx_RS_strong, idx_RS(i)];
end
end
RS_strong2 = unique(date(idx_RS_strong));


l5a_cell2 = unique(date(idx_l5a));
idx_l5a_strong = [];
for i = 1:length(idx_l5a)
if ~isempty(find(datapool.strong_idx==idx_l5a(i)))
idx_l5a_strong = [idx_l5a_strong, idx_l5a(i)];
end
end
l5a_strong2 = unique(date(idx_l5a_strong));
%% combine both for cells that ever generated bursts
IB_burst = unique([IB_burst1, IB_burst]);
RS_burst = unique([RS_burst1, RS_burst]);
l5a_burst = unique([l5a_burst1, l5a_burst]);
IB_cell = unique([IB_cell1, IB_cell]);
RS_cell = unique([RS_cell1, RS_cell]);
l5a_cell = unique([l5a_cell1, l5a_cell]);
l5a_cell{5} = '0910211';
figure
subplot(131), pie([length(l5a_cell)-length(l5a_burst),length(l5a_burst)])
subplot(132), pie([length(RS_cell)-length(RS_burst),length(RS_burst)])
subplot(133), pie([length(IB_cell)-length(IB_burst),length(IB_burst)])
% save(fullfile(pwd, 'stat', 'FigureS6_burst.mat'), 'l5a_cell', 'l5a_burst', 'RS_cell', 'RS_burst', 'IB_cell', 'IB_burst');
%% combine both for cells that ever generated spikelets
IB_strong = unique([IB_strong2, IB_strong1, IB_strong]);
RS_strong = unique([RS_strong2, RS_strong1, RS_strong]);
l5a_strong = unique([l5a_strong2, l5a_strong1, l5a_strong]);
IB_cell = unique([IB_cell2, IB_cell1, IB_cell]);
RS_cell = unique([RS_cell2, RS_cell1, RS_cell]);
l5a_cell = unique([l5a_cell2, l5a_cell1, l5a_cell]);
l5a_cell{5} = '0910211';
l5a_strong{3} = '0910211';
figure
subplot(131), pie([length(l5a_cell)-length(l5a_strong),length(l5a_strong)])
subplot(132), pie([length(RS_cell)-length(RS_strong),length(RS_strong)])
subplot(133), pie([length(IB_cell)-length(IB_strong),length(IB_strong)])

% save(fullfile(pwd, 'stat', 'FigureS6_Na.mat'), 'l5a_cell', 'l5a_strong', 'RS_cell', 'RS_strong', 'IB_cell', 'IB_strong');