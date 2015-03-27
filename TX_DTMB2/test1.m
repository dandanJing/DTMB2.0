PN_total_len = 432;
channelFilter = multipath_new(16,1/7.56,1,0);
channel_real = zeros(1,PN_total_len);
channel_real(1:length(channelFilter)) = channelFilter;
chan_freq = fft(channel_real,3888*8);
chan_temp = chan_freq(96:96:end);
len = length(chan_temp);
chan_ifft_temp = zeros(1,PN_total_len);
chan_ifft_temp(1:len) = ifft(chan_temp).*exp(j*2*pi/(3888*8)*(96-1)*(0:len-1));
chan_off = chan_ifft_temp./channel_real;
mse = norm(chan_ifft_temp-channel_real)/norm(channel_real)
figure;
plot(abs( chan_ifft_temp));
