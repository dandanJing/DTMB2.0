% PN_total_len = 432;
% channelFilter = multipath_new(16,1/7.56,1,0);
%     channel_real = zeros(1,PN_total_len);
%     channel_real(1:length(channelFilter)) = channelFilter;
%     chan_freq = fft(channel_real,2048);
%     for kk = 1:8
%         chan_temp = chan_freq(kk:8:2048);
%         chan_ifft_temp = zeros(1,PN_total_len);
%         chan_ifft_temp(1:256) = ifft(chan_temp).*exp(j*2*pi/2048*(kk-1)*(0:255));
%         chan_off = chan_ifft_temp./channel_real;
%         mse(kk) = norm(chan_ifft_temp-channel_real)/norm(channel_real);
%     end
%     figure;
%     plot(abs( chan_ifft_temp));

for kk = 50:100
    temp(kk,:) = sum(channel_estimate_spn(kk:kk+5,:));
end
figure;
plot(abs(temp(51,:)));