function chlen =chan_len_estimate_1(h_in)

abs_h_in = abs(h_in);
max_h = max(abs_h_in);
thresh = max_h * 0.1;
k = max(find(abs_h_in >= thresh));
chlen = k;
