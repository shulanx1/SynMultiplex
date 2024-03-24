close all
% clear all

clear V_sub Istep_sub out Vm Im length
idx_f = 2;
path = 'E:\data\DNAnanopore\022423\cell2_positive';
% path = 'E:\data\uncaging\physiological characterization\112120\cell5';
filename = sprintf('220226_001.sub and supra.%d.wcp', idx_f);
out=import_wcp(fullfile(path, filename),'debug');

clear length
n_channel = out.channel_no;
n_recording = length(out.rec_index);

dt = out.T(2)-out.T(1); % sample time in s
t = (-500*dt:dt:2500*dt)*1e3;

spike_time = {}; % spike time sorted as the rising phase when cross threshold
thr = 0; % voltage treshold for spike detection, in mV
Istep_sub = -50:10:10*(n_recording-1)-50;
Istep_supra = -200:50:50*(n_recording-1)-200;
V_reststate = [];
for i = 1:n_recording
    Vm(:,i) = out.S{3}(:,i);
    Im(:,i) = out.S{4}(:,i);
    spike_time{i} = [];
    for n = floor(2/dt)+1:floor(2.5/dt)
        if Vm(n, i)>thr && Vm(n-1, i)<=thr
            spike_time{i} = [spike_time{i}, n*dt]; % spike time in s
        end
    end
    Ihold = mean(Im(1:floor(0.5/dt),i)); % holding current, in pA
%     Istep_sub(i) =(mean(Im(floor(0.5/dt)+1:floor(1/dt),i)) - Ihold); %floor((mean(Im(floor(0.5/dt)+1:floor(1/dt),i)) - Ihold)/10)*10; % injected step current, in pA
%     Istep_supra(i) = floor((mean(Im(floor(2/dt)+1:floor(2.5/dt),i)) - Ihold)/50)*50; % injected step current, in pA
    V_sub(i) = mean(Vm(floor(0.7/dt)+1:floor(1/dt),i));
    Vrest(i) = mean(Vm(1:floor(0.5/dt) ,i));
    V_reststate = [V_reststate;Vm(1:floor(0.5/dt) ,i)];
%     if Istep_sub(i)== 0
%         Rin(i) = 0;
%     else
%         Rin(i) = (V_sub(i)-Vrest(i))/Istep_sub(i)*1000; % input resistance, MOhm
%     end
    FR(i) = length(spike_time{i})/0.5; % firing rate, Hz
end

Ihold = sum(mean(Im(1:floor(0.5/dt),:))/size(Im, 2)); % holding current, in pA
Vrest = sum(mean(Vm(1:floor(0.5/dt),:))/size(Vm,2));
% Rin_pooled = Rin(1);%mean(Rin(find(Istep_sub<-10))); % input resistance, MOhm
% Rin_pooled = (mean(Vm(floor(2.2/dt):floor(2.5/dt),1))-Vrest(1))/(Istep_supra(1))*1e3;
P = polyfit(Istep_sub([2:6]), V_sub([2:6]),1)*1e3;   % input resistance, MOhm
Rin_pooled = P(1);
 figure
scatter(V_sub, Istep_sub)

length = size(Im,2);
red = [0,0,0];
pink = [0.8,0.8,0.8];%[255, 192, 203]/255;
colors_p = [linspace(pink(1),red(1),length)', linspace(pink(2),red(2),length)', linspace(pink(3),red(3),length)'];

t = (-0.5:dt:1)*1e3;
% figure
% hold on
% for i = 1:size(Im,2)
%     plot(t, Im(floor(1.5/dt):floor(3/dt),i),'Color', colors_p(i,:))
% end

figure
hold on
for i = 1:size(Im,2)
plot(t, Vm(floor(1.5/dt):floor(3/dt),i),'Color', colors_p(i,:))
end
xlabel('t (ms)')
ylabel('V (mV)')

Rin_pooled
rms_vrest = rms(V_reststate-Vrest)
Vrest
save(fullfile(path,sprintf('sub_and_supra_%d.mat', idx_f)), 'Rin_pooled', 'Vm', 'Im', 'Istep_supra', 'FR', 't','Ihold','Vrest')
%%
i = 1;
data(i).date = '121922';
data(i).cell = '3';
data(i).El = Vrest;
data(i).Ee = 0;
data(i).Ei = -80;
%% pool
Im = [-200:50:250];
idx = [2,3];
for i = 1:3
FR_pool_mat{i} = cell2mat(FR_pool(idx,i));
end
figure
errorbar(Im,mean(FR_pool_mat{1}), std(FR_pool_mat{1})/sqrt(size(FR_pool_mat{1},1)), 'o', 'MarkerSize',10,'Color'  ,'k','MarkerEdgeColor', 'k','MarkerFaceColor', [1,1,1] )
hold on
errorbar(Im,mean(FR_pool_mat{2}), std(FR_pool_mat{2})/sqrt(size(FR_pool_mat{2},1)), 'o', 'MarkerSize',10,'Color',[1,0.5,0.5],  'MarkerEdgeColor', [1,0.5,0.5],'MarkerFaceColor', [1,1,1] )
errorbar(Im,mean(FR_pool_mat{3}), std(FR_pool_mat{3})/sqrt(size(FR_pool_mat{3},1)), 'o', 'MarkerSize',10,'Color', 'r', 'MarkerEdgeColor', 'r','MarkerFaceColor', [1,1,1] )
hold on
plot(Im,mean(FR_pool_mat{1}), 'Color', 'k')
plot(Im,mean(FR_pool_mat{2}),'Color',[1,0.5,0.5])
plot(Im,mean(FR_pool_mat{3}), 'Color', 'r')
ylabel('f [Hz]')
xlabel('I [pA]')
