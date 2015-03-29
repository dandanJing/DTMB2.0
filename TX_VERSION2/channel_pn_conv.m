function h_pn_conv_out = channel_pn_conv(PN,channel,chan_len)

PN_freq = fft(PN,2048);
PN_freq = reshape(PN_freq, 2048, 1);
chan_freq = fft(channel,2048);
chan_freq = reshape(chan_freq, 2048, 1);
chan_pn_conv = ifft(PN_freq.*chan_freq); 

total_len = length(PN)+chan_len;
h_pn_conv_out = zeros(1,2048);
h_pn_conv_out(1:total_len)=chan_pn_conv(1:total_len);
