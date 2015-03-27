function channel_smooth =chan_estimate_smooth(h_old,h_current,len,debug)

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
channel_smooth = h_temp/count_frame;