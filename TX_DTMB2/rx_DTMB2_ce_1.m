%%channel estimation 
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48*8, 64QAM
clear all,close all,clc
debug = 0;
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 8;%定义多径类型
SNR = [30];

%%参数定义
PN_total_len = 432; %帧头长度,前同步88，后同步89
DPN_total_len = 1024;
DPN_len = 512;
load pn256_pn512.mat
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
sim_num= 400; %仿真的帧数
iter_num = 1; %迭代次数
MAX_CHANNEL_LEN = PN_total_len;

%%帧头信号
PN = PN_gen*1.975;
temp = ifft(pn512);
DPN = temp*sqrt(var(PN)/var(temp));
dpn_h_smooth_alpha = 1/4;
dpn_h_smooth_result = [];
dpn_h_denoise_alpha = 1/256;

%%信道
if debug_multipath
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
    channel_real = zeros(1,PN_total_len);
    channel_real(1:length(channelFilter)) = channelFilter;
else
    channel_real = zeros(1,PN_total_len);
    channel_real(1) = 1;
end

dpn_start_frame = 30;
dpn_end_frame = sim_num-9;
spn_mean_mse = zeros(1,length(SNR));
dpn_mean_mse = zeros(1,length(SNR));
spn_pnrm_SNR = zeros(1,length(SNR));
dpn_pnrm_SNR = zeros(1,length(SNR));
spn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));
dpn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));

dpn_channel_mse = zeros(length(SNR),sim_num);
dpn_pn_rm_snr = zeros(length(SNR),sim_num);
dpn_chan_conv_snr = zeros(length(SNR),sim_num);
mse_pos = 1;
for SNR_IN = SNR %定义输入信噪比
    %%数据输入
    close all;
    figure;
    plot(abs(channel_real));
    title('真实信道估计');
    
    if debug_multipath
        matfilename = strcat('DTMB_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    else
        load strcat('DTMB_data_awgn_SNR',num2str(SNR_IN),'.mat');
    end
    
     %% 双PN信道估计
      %真实数据与真实信道的卷积结果
     channel_real_freq = fft(channel_real,FFT_len);
     frame_test = DPN_total_len + FFT_len;
     coeff = 6.4779e+04;
     for i=1:sim_num-1
          Receive_data = Send_data_srrc_tx1_dpn((i-1)*frame_test+(1:frame_test));
          pn_test = Receive_data(DPN_len+(1:DPN_len));
          pn_test = pn_test ./ coeff;
          pn512_fft = fft(pn_test);
          dpn_h_freq =  pn512_fft./ pn512;
          dpn_h_time = ifft(dpn_h_freq);
          if debug
              figure;
              plot(abs(dpn_h_time(1:PN_total_len)))
              title('双PN估计结果');
              pause;
          end
          chan_len_dpn = min(chan_len_estimate(dpn_h_time),MAX_CHANNEL_LEN);
          dpn_h_time(chan_len_dpn+1:end)=0;
          channel_estimate_dpn(i,1:PN_total_len)=dpn_h_time(1:PN_total_len);
          if i==1
              dpn_h_smooth_result = dpn_h_time(1:PN_total_len);
          else
              dpn_h_current_frame = mean(channel_estimate_dpn(i-1:i,1:PN_total_len));
              dpn_h_smooth_result = dpn_h_smooth_alpha*dpn_h_time(1:PN_total_len)+(1-dpn_h_smooth_alpha)*dpn_h_smooth_result;
          end
          %dpn_h_smooth_result = channel_denoise2(dpn_h_smooth_result,dpn_h_denoise_alpha);
          chan_len_dpn = min(chan_len_estimate(dpn_h_smooth_result),MAX_CHANNEL_LEN);
          dpn_h_smooth_result(chan_len_dpn+1:end)=0;
          dpn_channel_mse(mse_pos,i) = norm(dpn_h_smooth_result-channel_real)/norm(channel_real);

          %%数据恢复
          dpn_h_freq_frame = fft(dpn_h_smooth_result,FFT_len);

          pn_conv = channel_pn_conv( DPN, dpn_h_smooth_result,chan_len_dpn);
          frame_tail =  Send_data_srrc_tx1_dpn(i*frame_test+(1:chan_len_dpn))-pn_conv(1:chan_len_dpn);
          frame_data =  Receive_data(DPN_total_len+1:end);
          frame_data(1:chan_len_dpn) = frame_data(1:chan_len_dpn)-pn_conv(DPN_len+(1:chan_len_dpn))+frame_tail;

          %真实信道与真实数据的卷积后的时域和频域结果
          data_real_freq =  data_transfer((i-1)*FFT_len+(1:FFT_len)).*channel_real_freq;
          data_real_time = ifft(data_real_freq);
          
          dpn_pn_rm_snr(mse_pos,i) = estimate_SNR(frame_data,data_real_time);
          dpn_temp = data_transfer((i-1)*FFT_len+(1:FFT_len)).*dpn_h_freq_frame;
          dpn_chan_conv_snr(mse_pos,i) = estimate_SNR(dpn_temp,data_real_freq);
     end
     fft_data_dpn = fft(frame_data)./channel_real_freq;
      figure;
      plot(fft_data_dpn,'.k');
      title('双PN数据均衡结果');
      
      dpn_mean_mse(mse_pos) = mean(dpn_channel_mse(mse_pos,dpn_start_frame:dpn_end_frame))
      dpn_pnrm_SNR(mse_pos) = mean(dpn_pn_rm_snr(mse_pos,dpn_start_frame:dpn_end_frame))
      dpn_pn_chan_conv_freq_SNR(mse_pos) = mean(dpn_chan_conv_snr(mse_pos,dpn_start_frame:dpn_end_frame))
      mse_pos = mse_pos +1;
end

figure;
plot(abs(dpn_h_smooth_result));
title('双PN信道估计结果');