clear all,close all,clc

debug = 0;
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 18;%定义多径类型
SNR = [25:5:35];
debug_TPS_mode = 2;

%%参数定义
PN_total_len = 432; %帧头长度,前同步88，后同步89
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
sim_num= 1000; %仿真的帧数
MAX_CHANNEL_LEN = PN_total_len;
load TPS_symbol.mat

%%真实信道
if debug_multipath
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
    channel_real = zeros(1,PN_total_len);
    channel_real(1:length(channelFilter)) = channelFilter;
else
    channel_real = zeros(1,PN_total_len);
    channel_real(1) = 1;
end
channel_real_freq = fft(channel_real,FFT_len);

spn_start_measure_mse_frame = 50;
spn_end_measure_mse_frame = sim_num-9;
spn_h_freq_off_thresh = 0.02;
spn_h_denoise_alpha = 0.0005;
spn_h_smooth_alpha = 1/4;

mse_pos = 1;
for SNR_IN = SNR %定义输入信噪比
    %%数据输入
    close all;
    
    if debug_multipath
        matfilename = strcat('DTMB2_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    else
        matfilename = strcat('DTMB2_data_awgn_SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    end
    
    %%单PN信道估计,干扰消除
    start_pos = 0;
    chan_len_spn = PN_total_len;
    last_frame_tail = zeros(1,chan_len_spn);
    frame_tps_div_scatter_delay = [];
    frame_tps_div_scatter_delay2 = [];
    last_frame_tps_ifft = [];
    spn_h_smooth_result = zeros(1,PN_total_len);
    channel_estimate_spn = zeros(sim_num,PN_total_len);
    for i=1:sim_num-1
          close all;
          if debug || i == sim_num-1
              figure;
              plot(abs(channel_real));
              title('真实信道估计');
          end   
          
          Receive_data = Send_data_srrc_tx1_spn(start_pos+(1:Frame_len));
          pn_frame_temp = Receive_data(1:PN_total_len);
          pn_frame_temp(1:chan_len_spn) = pn_frame_temp(1:chan_len_spn)-last_frame_tail;
          
          PN = PN_gen(i,0);
          h_pn_conv = channel_pn_conv(PN, spn_h_smooth_result, chan_len_spn);
          PN_next = PN_gen(i+1,0);
          h_pn_conv_next = channel_pn_conv(PN_next, spn_h_smooth_result, chan_len_spn);
          frame_data_time =  Receive_data(PN_total_len+(1:FFT_len));
          pn_tail = h_pn_conv(PN_total_len+(1:chan_len_spn-1));
          data_tail =  Send_data_srrc_tx1_spn(start_pos+Frame_len+(1:chan_len_spn-1))-h_pn_conv_next(1:chan_len_spn-1);
          frame_data_recover =  frame_data_time;
          frame_data_recover(1:chan_len_spn-1)=frame_data_recover(1:chan_len_spn-1)-pn_tail+data_tail;
          frame_data_recover_freq = fft(frame_data_recover);  
          
          %%TPS channel estimate
          current_tps_position = TPS_pos_block81_gen(i,0);
          current_tps_div = frame_data_recover_freq(current_tps_position)./tps_value;
          
          if debug_TPS_mode == 1
              h_iter = ifft(current_tps_div);
              h_iter = h_iter.* exp(1j*2*pi/(FFT_len)*(0:length(h_iter)-1)*(current_tps_position(1)-1));
          elseif debug_TPS_mode == 2
              current_tps_ifft = ifft(current_tps_div);
              current_tps_ifft = current_tps_ifft.* exp(1j*2*pi/(FFT_len)*(0:length(current_tps_ifft)-1)*(current_tps_position(1)-1));
              last_tps_pos = TPS_pos_block81_gen(i-1,0);
              h_iter = zeros(1,PN_total_len);
              if i == 1
                  last_frame_tps_ifft = current_tps_ifft;
              end
              h_iter(1:384)=(last_frame_tps_ifft+current_tps_ifft)/2;
              h_iter(385:end)=(current_tps_ifft(1:432-384)-last_frame_tps_ifft(1:432-384))/...
                  (exp(-j*2*pi*384/3888/8*(current_tps_position(1)-1))-exp(-j*2*pi*384/3888/8*(last_tps_pos(1)-1)));
              h_iter(1:432-384)=h_iter(1:432-384)-h_iter(385:end)*(exp(-j*2*pi*384/3888/8*(current_tps_position(1)-1))+exp(-j*2*pi*384/3888/8*(last_tps_pos(1)-1)))/2;
          elseif debug_TPS_mode == 3
              tps_pos_total = [tps_pos_scatter TPS_pos_block81_gen(i-1,0) TPS_pos_block81_gen(i-2,0)];
              tps_value_total = [tps_scatter_div frame_tps_div_scatter_delay frame_tps_div_scatter_delay2];
              tps_temp(tps_pos_total)=tps_value_total;
              tps_pos_total = sort(tps_pos_total);
              tps_value_total = tps_temp(tps_pos_total);
              h_iter = ifft(tps_value_total);
              h_iter = h_iter.* exp(1j*2*pi/(FFT_len)*(0:length(h_iter)-1)*(tps_pos_total(1)-1));
          end
          
          h_iter = channel_denoise(h_iter, spn_h_denoise_alpha);
          h_temp = zeros(1,PN_total_len);
          len = min(PN_total_len,length(h_iter));
          h_temp(1:len)=h_iter(1:len);
          h_iter = h_temp;
          clear h_temp;
           
          %channel length
          if i <= 8
             [h_iter chan_len_spn] = chan_estimate(channel_estimate_spn(1:i-1,:),h_iter,i-1,0);
          else
             [h_iter chan_len_spn] = chan_estimate(channel_estimate_spn(i-8:i-1,:),h_iter,8,0);
             h_iter = spn_h_smooth_alpha*h_iter+(1-spn_h_smooth_alpha)*spn_h_smooth_result;
          end
          
           % save status
           if debug_TPS_mode == 3
             frame_tps_div_scatter_delay2 = frame_tps_div_scatter_delay;
             frame_tps_div_scatter_delay = current_tps_div;
           elseif debug_TPS_mode == 2
               last_frame_tps_ifft = current_tps_ifft;
           end
          
          channel_len_estimate_spn(mse_pos,i) = chan_len_spn;
          channel_estimate_spn(i,:) = h_iter;
          spn_h_smooth_result = h_iter;
          
          if debug
              figure;
              plot(abs(h_iter));
              title('TPS估计结果');
              pause;
          end
          %%数据干扰消除
          h_pn_conv = channel_pn_conv(PN, h_iter, chan_len_spn);
          h_pn_conv_next = channel_pn_conv(PN_next, h_iter, chan_len_spn);
          spn_h_freq = fft(h_iter,FFT_len);                 
          frame_data_time =  Receive_data(PN_total_len+(1:FFT_len));
          pn_tail = h_pn_conv(PN_total_len+(1:chan_len_spn-1));
          data_tail =  Send_data_srrc_tx1_spn(start_pos+Frame_len+(1:chan_len_spn-1))-h_pn_conv_next(1:chan_len_spn-1);
          frame_data_recover =  frame_data_time;
          frame_data_recover(1:chan_len_spn-1)=frame_data_recover(1:chan_len_spn-1)-pn_tail+data_tail;
          
          %信道估计 mse
          spn_channel_mse(mse_pos,i) = norm(spn_h_smooth_result-channel_real)/norm(channel_real);
          
          %计算当前帧的拖尾，以在下一帧消除其对于PN估计的影响
          frame_data_eq_freq = fft(frame_data_recover)./spn_h_freq;
          frame_data_eq_time = ifft(frame_data_eq_freq);
          frame_data_tail =  frame_data_eq_time(FFT_len-PN_total_len+(1:PN_total_len));
          data_tail_chan_conv = channel_pn_conv(frame_data_tail,spn_h_smooth_result, chan_len_spn);
          last_frame_tail = data_tail_chan_conv(PN_total_len+(1:chan_len_spn));
          
          %真实信道与真实数据的卷积后的时域和频域结果
          data_real_freq =  data_transfer((i-1)*FFT_len+(1:FFT_len)).*channel_real_freq;
          data_real_time = ifft(data_real_freq);

          spn_pn_rm_snr(mse_pos,i) = estimate_SNR(frame_data_recover,data_real_time);
          spn_temp = data_transfer((i-1)*FFT_len+(1:FFT_len)).*spn_h_freq;
          spn_chan_conv_snr(mse_pos,i) = estimate_SNR(spn_temp,data_real_freq);

          start_pos = start_pos + Frame_len;
    end
    max_mse(mse_pos) = max(spn_channel_mse(mse_pos,spn_start_measure_mse_frame:spn_end_measure_mse_frame));
    min_mse(mse_pos) = min(spn_channel_mse(mse_pos,spn_start_measure_mse_frame:spn_end_measure_mse_frame));
    spn_mean_mse(mse_pos) = mean(spn_channel_mse(mse_pos,spn_start_measure_mse_frame:spn_end_measure_mse_frame));
    spn_pnrm_SNR(mse_pos) = mean(spn_pn_rm_snr(mse_pos,spn_start_measure_mse_frame:spn_end_measure_mse_frame));
    spn_pn_chan_conv_freq_SNR(mse_pos) = mean(spn_chan_conv_snr(mse_pos,spn_start_measure_mse_frame:spn_end_measure_mse_frame));
    temp1 = 1/(10^(spn_pnrm_SNR(mse_pos)/10))+ 1/(10^(spn_pn_chan_conv_freq_SNR(mse_pos)/10));
    spn_ce_SNR(mse_pos) = 10*log10(1/temp1);
    mse_pos = mse_pos +1;
end
max_mse
min_mse
spn_mean_mse
spn_pnrm_SNR
spn_pn_chan_conv_freq_SNR
spn_ce_SNR

figure;
plot(abs(spn_h_smooth_result));
title('单PN信道估计结果');

figure;
semilogy(SNR,spn_mean_mse,'r*-');
title('单PN估计MSE');

figure;
plot(SNR,spn_pnrm_SNR,'r*-');
title('数据循环重构后的信噪比');

figure;
plot(SNR,spn_pn_chan_conv_freq_SNR,'r*-');
title('信道估计误差的信噪比');