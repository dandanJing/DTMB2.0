function chlen =chan_len_estimate_new(h_in,alpha)

abs_h_in = abs(h_in);
max_h = max(abs_h_in);
thresh = max_h * alpha;
k = max(find(abs_h_in >= thresh));
chlen = k+20;
