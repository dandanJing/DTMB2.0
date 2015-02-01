function SNR_out = estimate_SNR(estimate, data_source)
len = length(data_source);
noise = (data_source-estimate);
noise_ave = mean(noise);
noise_power = sum((abs(noise-noise_ave)).^2)/len;
source_ave = mean(data_source);
source_power = sum(abs(data_source-source_ave).^2)/len;
SNR_out =10*log( source_power / noise_power);