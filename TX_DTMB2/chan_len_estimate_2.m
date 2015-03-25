function chlen =chan_len_estimate_2(h_old,h_current,len)

count_frame = 1;
h_temp = h_current;
len_temp = len;
if len > 8
    len_temp = 8;
end
for kk = 1:len_temp
    h_old_temp = h_old(kk,:);
    if sum(h_current.*conj(h_old_temp))/sum(h_current.*conj(h_current)) > 0.9
        count_frame = count_frame + 1;
        h_temp = h_temp + h_old_temp;
    end
end
abs_h_temp = abs(h_temp);
max_h = max(abs_h_temp);
alpha = 0.33 - (count_frame-1)*0.04;
thresh = max_h * alpha;
chlen = max(find(abs_h_temp >= thresh));
