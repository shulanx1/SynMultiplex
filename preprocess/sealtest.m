%% Ch2
close all
% clear all

idx = 4;  
path = 'E:\data\DNAnanopore\030223\cell4_acsf_qx314';
filename = sprintf('230303_001.sealtest.%d.wcp', idx(1));
out=import_wcp(fullfile(path, filename),'debug');

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{3};
Vm = out.S{4};
t = 0:dt:dt*(size(Vm,1)-1);

I_baseline = mean(Im(1:4000,1));
I_step = mean(Im(6000:12000,1))-I_baseline;
I_trans = I_baseline-min(Im(:,1));
Cm = sum((Im(4000:6000,1)-I_baseline)-I_step)*dt*1e3/10; %pF

Ra = 10/(I_trans/1e3);  %in Mohm
Rm = 10/(I_step/1e3)-Ra;
Rseal = -10/(I_step/1e3)
figure
plot(t, Im(:,3))
%% Ch1
close all
% clear all

idx = 7;  
path = 'E:\data\DNAnanopore\011324\cell6_wo';
filename = sprintf('240114_001.sealtest_Ch1.%d.wcp', idx(1));
out=import_wcp(fullfile(path, filename),'debug');

n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s

   
Im = out.S{1};
Vm = out.S{2};
t = 0:dt:dt*(size(Vm,1)-1);

I_baseline = mean(Im(1:4000,1));
I_step = mean(Im(6000:12000,1))-I_baseline;
I_trans = I_baseline-min(Im(:,1));
Cm = sum((Im(4000:6000,1)-I_baseline)-I_step)*dt*1e3/10; %pF

Ra = 10/(I_trans/1e3);  %in Mohm
Rm = -10/(I_step/1e3)-Ra;
Rseal = -10/(I_step/1e3)
figure
plot(t, Im(:,1))
%%
i = 1;
data(i).date = '121922';
data(i).cell = '3';
data(i).Cm = Cm; %pF
data(i).gl = 1e3/Rm; %nS
data(i).Ee = 0;
data(i).Ei = -80;
%%
figure
plot(t, Im(:,1))

