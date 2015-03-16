function h_out = channel_filter(channel_in,thresh)
h_out = channel_in;
h_max = max(abs(h_out));
h_temp = channel_in;
len = length(channel_in);
ave_len = 5;
pos_all = find(abs(h_out)<h_max*thresh);
for kk = pos_all
    pos_temp = find(pos_all>kk-ave_len);
    pos1 = pos_temp(1);
    pos_temp = find(pos_all<kk+ave_len);
    pos2 = pos_temp(end);
    h_out(kk) = mean(h_temp(pos_all(pos1:pos2)));
end
% h_max = max(abs(h_out));
% h_out(abs(h_out)<h_max*thresh)=h_out(abs(h_out)<h_max*thresh)*0.5;