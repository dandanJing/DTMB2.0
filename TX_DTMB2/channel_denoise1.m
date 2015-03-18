function h_out = channel_denoise1(channel_in,thresh)

channel_mmse = channel_mmse_filter(channel_in,0.99).';
h_out = channel_mmse;
h_max = max(abs(channel_in));
h_out(abs(channel_in)>h_max*thresh)=channel_in(abs(channel_in)>h_max*thresh);
