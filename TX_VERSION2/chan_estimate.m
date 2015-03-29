function [channel chlen] =chan_estimate(h_old,h_current,len,debug)

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
alpha = 0.33 - (count_frame-1)*0.03;
thresh = max_h * alpha;
chlen = max(find(abs_h_temp >= thresh));
if chlen < length(h_current)
    thresh_num = length(find(abs(h_current)>0.15*max(abs(h_current))));
    if thresh_num > 1
        chlen = chlen + 20;
        chlen = min(chlen,length(h_current));
    end
end 
channel = h_temp/count_frame;
channel(chlen+1:end) = 0;
if debug
    chlen
    alpha
    count_frame
    figure;
    plot(abs(abs_h_temp));
    title('PN长度估计时的累加结果');
end