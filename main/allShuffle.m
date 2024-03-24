function Total_shuffled_spiketime = allShuffle(SpikeTime, t,surrogate)
if nargin<3 || strcmp(surrogate,''), surrogate = 5; end
Spikes = zeros(length(SpikeTime), length(t));
for n = 1:length(SpikeTime)
    Spikes(n,:) = ismember(t,SpikeTime{n});
end
cell_count = size(Spikes,1);
num_images = size(Spikes,2);
Total_shuffled = zeros(cell_count,num_images);
for x = 1:surrogate
    p = randperm(num_images);
    q = randperm(cell_count)';
    for i = 1:cell_count
        for j = 1:num_images
            Total_shuffled(i,j) = Spikes(q(i),p(j));
        end
    end
end
Total_shuffled_spiketime = [];
for n = 1:length(SpikeTime)
    Total_shuffled_spiketime = [Total_shuffled_spiketime,t(find(Total_shuffled(n,:)~=0))];
end
end

