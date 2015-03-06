function h_out = channel_denoise(channel_in,thresh)
h_out = channel_in;
h_max = max(abs(h_out));
h_out(abs(h_out)<h_max*thresh)=h_out(abs(h_out)<h_max*thresh)*0.5;
