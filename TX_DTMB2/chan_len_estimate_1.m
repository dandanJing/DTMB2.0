function chlen =chan_len_estimate_1(h_temp)

h_in = channel_mmse_filter(h_temp, 0.99);
abs_h_in = abs(h_in);
max_h = max(abs_h_in);
thresh = max_h * 0.05;
k = max(find(abs_h_in >= thresh));
chlen = k+20;
