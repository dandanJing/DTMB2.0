clear all,close all,clc
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 16;%定义多径类型
PN_total_len = 432;
SNR = [100];
debug_frame_two = 1;

%%真实信道
if debug_multipath
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
    channel_real = zeros(1,PN_total_len);
    channel_real(1:length(channelFilter)) = channelFilter;
else
    channel_real = zeros(1,PN_total_len);
    channel_real(1) = 1;
end

channel_real_awgn =  awgn(channel_real,SNR,'measured');
channel_real_SNR = 10*log10(var(channel_real)/var(channel_real_awgn-channel_real))
noise = channel_real_awgn-channel_real;
channel_freq = fft(channel_real,3888*8);
channel_freq_awgn = fft(channel_real_awgn,3888*8);
channel_freq_SNR = 10*log10(var(channel_freq)/var(channel_freq_awgn-channel_freq))
sim_num = 1;
modify_len = 432-384;
for mm = 1:sim_num
    for kk = 1:80
        tps_pos = [kk:81:3888*8];
        tps_value = channel_freq_awgn(tps_pos);
        chan_es = zeros(1,PN_total_len);
        len = min(PN_total_len,length(tps_value));
        chan_es_temp = ifft(tps_value);
        chan_es_temp = chan_es_temp .*exp(j*2*pi*(kk-1)*(0:length(tps_pos)-1)/(3888*8));
        chan_es(1:len) = chan_es_temp(1:len);
        if debug_frame_two
            if kk >1
                h_recover = zeros(1,PN_total_len);
                h_recover(1:len)=chan_es_temp(1:len);
                h_recover(48*8+1:end)=(chan_es(1:modify_len)-chan_es_old(1:modify_len))/(exp(-j*2*pi*(tps_pos(1)-1)/(3888*8)*384)-exp(-j*2*pi*(tps_pos_old(1)-1)/(3888*8)*384));
                h_recover(1:modify_len)= h_recover(1:modify_len)-h_recover(48*8+1:end)*exp(-j*2*pi*(tps_pos(1)-1)/(3888*8)*384);
                ce_mmse(mm,kk) =  norm(h_recover-channel_real)/norm(channel_real);
                SNR(mm,kk) = 10*log10(var(channel_real)/var(h_recover-channel_real));
            end
        else
            ce_mmse(mm,kk) = norm(chan_es-channel_real)/norm(channel_real);
            SNR(mm,kk) = 10*log10(var(channel_real)/var(chan_es-channel_real)); 
        end     
        tps_pos_old = tps_pos;
        chan_es_old = chan_es;       
    end
end

ce_mean = mean(ce_mmse)
SNR_mean = mean(SNR)
[ce_Y ce_I] = min(ce_mean);
[SNR_Y SNR_I]=max(SNR_mean);

figure;
plot(abs(channel_real));
title('真实信道估计');
if debug_frame_two
    figure;
    plot(abs(h_recover));
    title('信道估计结果');
else
    figure;
    plot(abs(chan_es));
    title('信道估计结果');
end

