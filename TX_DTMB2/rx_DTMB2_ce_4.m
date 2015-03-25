%%channel estimation 数据干扰消除 TPS
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48*8, 64QAM
clear all,close all,clc

debug = 0;
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 8;%定义多径类型
SNR = [25:5:35];

%%参数定义
PN_total_len = 432; %帧头长度,前同步88，后同步89
DPN_total_len = 1024;
DPN_len = 512;
load pn256_pn512.mat
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
sim_num= 1000; %仿真的帧数
iter_num = 2; %迭代次数
MAX_CHANNEL_LEN = PN_total_len;

%%帧头信号
PN = PN_gen*1.975;
PN_binary = conj(PN_gen)/1024;
temp = ifft(pn512);
DPN = temp*sqrt(var(PN)/var(temp));
dpn_h_smooth_alpha = 1/4;
dpn_h_smooth_result = [];
dpn_h_denoise_alpha = 1/256;
spn_h_denoise_alpha = 0.01;
spn_h_denoise_alpha_iter1 = 0.2;
spn_tps_h_denoise_alpha = 1/256;
spn_h_smooth_alpha = 1/4;
spn_h_smooth_result = zeros(1,PN_total_len);
spn_h_pn_smooth_result = zeros(1,PN_total_len);
spn_h_off_thresh = 0.02;
spn_h_len_alpha_pn = 0.03;
spn_h_len_alpha_TPS = 0.03;
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
dpn_start_frame = 50;
dpn_end_frame = sim_num-9;
spn_start_frame = 50;
spn_denoise_decrease_frame = 30;
spn_end_frame = sim_num-9;
spn_mean_mse = zeros(1,length(SNR));
dpn_mean_mse = zeros(1,length(SNR));
spn_pnrm_SNR = zeros(1,length(SNR));
dpn_pnrm_SNR = zeros(1,length(SNR));
spn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));
dpn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));
spn_ce_SNR = zeros(1,length(SNR));
dpn_ce_SNR = zeros(1,length(SNR));

dpn_channel_mse = zeros(length(SNR),sim_num);
dpn_pn_rm_snr = zeros(length(SNR),sim_num);
dpn_chan_conv_snr = zeros(length(SNR),sim_num);
spn_channel_mse = zeros(length(SNR),sim_num);
spn_pn_rm_snr = zeros(length(SNR),sim_num);
spn_chan_conv_snr = zeros(length(SNR),sim_num);

channel_estimate_dpn = zeros(length(SNR),PN_total_len);
channel_estimate_spn = zeros(length(SNR),PN_total_len);
mse_pos = 1;
for SNR_IN = SNR %定义输入信噪比
    %%数据输入
    close all;
    
    if debug_multipath
        matfilename = strcat('DTMB_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    else
        matfilename = strcat('DTMB_data_awgn_SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    end
    %Send_data_srrc_tx1_ch_spn 单PN过多径信道
    %Send_data_srrc_tx1_spn 单PN过多径和高斯信道
    %Send_data_srrc_tx1_ch_dpn 双PN过多径信道
    %Send_data_srrc_tx1_dpn 双PN过多径和高斯信道
    %data_transfer  传输数据符号
    %tps_position  TPS位置
    %tps_symbol TPS符号
    %tps_block_len 每帧中TPS块长度
    %dimY  几个符号构成一个完整的导频帧
    %dimX_len OFDM符号间TPS位置间隔
    %Super_Frame  超帧长度
    
     %% 双PN信道估计
      %真实数据与真实信道的卷积结果
     frame_test = DPN_total_len + FFT_len;
     coeff = 6.4779e+04;
     for i=1:sim_num-1
%           Receive_data = Send_data_srrc_tx1_dpn((i-1)*frame_test+(1:frame_test));
%           pn_test = Receive_data(DPN_len+(1:DPN_len));
%           pn_test = pn_test ./ coeff;
%           pn512_fft = fft(pn_test);
%           dpn_h_freq =  pn512_fft./ pn512;
%           dpn_h_time = ifft(dpn_h_freq);
%           dpn_h_time = channel_denoise2(dpn_h_time,dpn_h_denoise_alpha);
%           if debug
%               figure;
%               plot(abs(dpn_h_time(1:PN_total_len)))
%               title('双PN估计结果');
%               pause;
%           end
%           chan_len_dpn = min(chan_len_estimate(dpn_h_time),MAX_CHANNEL_LEN);
%           dpn_h_time(chan_len_dpn+1:end)=0;
%           channel_estimate_dpn(i,1:PN_total_len)=dpn_h_time(1:PN_total_len);
%           if i==1
%               dpn_h_smooth_result = dpn_h_time(1:PN_total_len);
%           else
%               dpn_h_current_frame = mean(channel_estimate_dpn(i-1:i,1:PN_total_len));
%               dpn_h_smooth_result = dpn_h_smooth_alpha*dpn_h_time(1:PN_total_len)+(1-dpn_h_smooth_alpha)*dpn_h_smooth_result;
%           end
%           %
%           chan_len_dpn = min(chan_len_estimate(dpn_h_smooth_result),MAX_CHANNEL_LEN);
%           dpn_h_smooth_result(chan_len_dpn+1:end)=0;
%           dpn_channel_mse(mse_pos,i) = norm(dpn_h_smooth_result-channel_real)/norm(channel_real);
% 
%           %%数据恢复
%           dpn_h_freq_frame = fft(dpn_h_smooth_result,FFT_len);
% 
%           pn_conv = channel_pn_conv( DPN, dpn_h_smooth_result,chan_len_dpn);
%           frame_tail =  Send_data_srrc_tx1_dpn(i*frame_test+(1:chan_len_dpn))-pn_conv(1:chan_len_dpn);
%           frame_data =  Receive_data(DPN_total_len+1:end);
%           frame_data(1:chan_len_dpn) = frame_data(1:chan_len_dpn)-pn_conv(DPN_len+(1:chan_len_dpn))+frame_tail;
% 
%           %真实信道与真实数据的卷积后的时域和频域结果
%           data_real_freq =  data_transfer((i-1)*FFT_len+(1:FFT_len)).*channel_real_freq;
%           data_real_time = ifft(data_real_freq);
%           
%           dpn_pn_rm_snr(mse_pos,i) = estimate_SNR(frame_data,data_real_time);
%           dpn_temp = data_transfer((i-1)*FFT_len+(1:FFT_len)).*dpn_h_freq_frame;
%           dpn_chan_conv_snr(mse_pos,i) = estimate_SNR(dpn_temp,data_real_freq);
     end
%       fft_data_dpn = fft(frame_data)./channel_real_freq;
%       if debug
%            figure;
%            plot(fft_data_dpn,'.k');
%            title('双PN数据均衡结果');
%       end
%      
%       dpn_mean_mse(mse_pos) = mean(dpn_channel_mse(mse_pos,dpn_start_frame:dpn_end_frame))
%       dpn_pnrm_SNR(mse_pos) = mean(dpn_pn_rm_snr(mse_pos,dpn_start_frame:dpn_end_frame))
%       dpn_pn_chan_conv_freq_SNR(mse_pos) = mean(dpn_chan_conv_snr(mse_pos,dpn_start_frame:dpn_end_frame))
%       temp2 = 1/(10^(dpn_pnrm_SNR(mse_pos)/10))+ 1/(10^(dpn_pn_chan_conv_freq_SNR(mse_pos)/10));
%       dpn_ce_SNR(mse_pos) = 10*log10(1/temp2);
      
      %%单PN信道估计,PN拼接
      start_pos = 0;
      h_iter = zeros(1,PN_total_len);
      chan_len_spn = PN_total_len;
      chan_len_spn_pn = PN_total_len;
      h_pn_conv = [];
      last_frame_tail = zeros(1,chan_len_spn);
      last_frame_h_tps = [];
      last2_frame_h_tps = [];
      for i=1:sim_num-1
          close all;
          if debug || i == sim_num-1
              figure;
              plot(abs(channel_real));
              title('真实信道估计');
          end   
          
          Receive_data = Send_data_srrc_tx1_spn(start_pos+(1:Frame_len));
          %PN
          pn_data = Receive_data(177+1:PN_total_len);
          pn_binary_temp = PN_binary(177+1:end);
          h_pn_es = zeros(1,255);
          for mm = 1:255
              h_pn_es(mm) = (sum(pn_data.*[pn_binary_temp(255-mm+2:end) pn_binary_temp(1:255-mm+1)]))/(1024*1.975*2*255);
          end
           h_pn_es = channel_denoise2(h_pn_es, spn_tps_h_denoise_alpha);
           if i ~= 1
                h_pn_es = spn_h_smooth_alpha*h_pn_es+(1-spn_h_smooth_alpha)*spn_h_pn_smooth_result;
           end
           if debug || i == sim_num-1
              figure;
              plot(abs(h_pn_es));
              title('PN原始估计结果');
          end
           chan_len_spn_pn_old =  chan_len_spn_pn;
           chan_len_spn_pn = min(chan_len_estimate_new(h_pn_es, spn_h_len_alpha_pn),MAX_CHANNEL_LEN);
           h_pn_es(chan_len_spn_pn+1:end)=0;
           spn_h_pn_smooth_result = h_pn_es;
           
          if debug || i == sim_num-1
              spn_h_len_alpha_pn
              chan_len_spn_pn
              figure;
              plot(abs(h_pn_es));
              title('PN估计结果');
          end
          
          %TPS
          h_pn_conv = channel_pn_conv(PN, spn_h_smooth_result, chan_len_spn);                  
          frame_data_time =  Receive_data(PN_total_len+(1:FFT_len));
          pn_tail = h_pn_conv(PN_total_len+(1:chan_len_spn));
          data_tail =  Send_data_srrc_tx1_spn(start_pos+Frame_len+(1:chan_len_spn))-h_pn_conv(1:chan_len_spn);
          frame_data_recover =  frame_data_time;
          frame_data_recover(1:chan_len_spn)=frame_data_recover(1:chan_len_spn)-pn_tail+data_tail;
          frame_data_recover_freq = fft(frame_data_recover);
          tps_pos_current = tps_position+mod(i,dimY)*dimX_len; 
          h_frame_tps_current = frame_data_recover_freq(tps_pos_current)./tps_symbol;
          if i > 2
              tps_pos_total = [tps_pos_current tps_position+mod(i-1,dimY)*dimX_len tps_position+mod(i-2,dimY)*dimX_len];
             h_tps_total = [h_frame_tps_current last_frame_h_tps last2_frame_h_tps];
          else
              tps_pos_total =  tps_pos_current;
              h_tps_total =  h_frame_tps_current;
          end
          h_frame_tps(tps_pos_total) =  h_tps_total;
          pos_temp = sort(tps_pos_total);
          h_tps_modify = h_frame_tps(pos_temp);
          h_tps_es = ifft(h_tps_modify);
          h_tps_es_total = h_tps_es;
          len_temp = min(PN_total_len,length(h_tps_es));
          h_tps_es = zeros(1,PN_total_len);
          h_tps_es(1:len_temp) = h_tps_es_total(1:len_temp);
          h_tps_es = channel_denoise2(h_tps_es, spn_tps_h_denoise_alpha);
          if i~= 1
              h_tps_es = spn_h_smooth_alpha *h_tps_es+(1-spn_h_smooth_alpha)*spn_h_smooth_result;
          end
          chan_len_spn_old =  chan_len_spn;
          chan_len_spn = min(chan_len_estimate_new(h_tps_es, spn_h_len_alpha_TPS),MAX_CHANNEL_LEN);
           h_tps_es(chan_len_spn+1:end) = 0;
          spn_h_smooth_result = h_tps_es;
          if debug
              chan_len_spn
              figure;
              plot(abs(h_tps_es));
              title('单PN TPS估计结果');
              pause;
          end
          
          h_result = spn_h_smooth_result;
          if chan_len_spn < 200
              h_result(1:length(spn_h_pn_smooth_result)) = spn_h_pn_smooth_result;
          end  
          channel_estimate_spn(i,:)= h_result;
          last2_frame_h_tps = last_frame_h_tps;
          last_frame_h_tps = h_frame_tps_current;
          
          spn_h_freq = fft(h_result,FFT_len);
          spn_channel_mse(mse_pos,i) = norm(h_result-channel_real)/norm(channel_real);
          %计算当前帧的拖尾，以在下一帧消除其对于PN估计的影响
          frame_data_eq_freq = fft(frame_data_recover)./spn_h_freq;
          frame_data_eq_time = ifft(frame_data_eq_freq);
          frame_data_tail =  frame_data_eq_time(FFT_len-PN_total_len+(1:PN_total_len));
          data_tail_chan_conv = channel_pn_conv(frame_data_tail,h_result, chan_len_spn);
          last_frame_tail = data_tail_chan_conv(PN_total_len+(1:chan_len_spn));
          
          %真实信道与真实数据的卷积后的时域和频域结果
          data_real_freq =  data_transfer((i-1)*FFT_len+(1:FFT_len)).*channel_real_freq;
          data_real_time = ifft(data_real_freq);

          spn_pn_rm_snr(mse_pos,i) = estimate_SNR(frame_data_recover,data_real_time);
          spn_temp = data_transfer((i-1)*FFT_len+(1:FFT_len)).*spn_h_freq;
          spn_chan_conv_snr(mse_pos,i) = estimate_SNR(spn_temp,data_real_freq);

          start_pos = start_pos + Frame_len;
      end
      spn_mean_mse(mse_pos) = mean(spn_channel_mse(mse_pos,spn_start_frame:spn_end_frame))
      spn_pnrm_SNR(mse_pos) = mean(spn_pn_rm_snr(mse_pos,spn_start_frame:spn_end_frame))
      spn_pn_chan_conv_freq_SNR(mse_pos) = mean(spn_chan_conv_snr(mse_pos,spn_start_frame:spn_end_frame))
      temp1 = 1/(10^(spn_pnrm_SNR(mse_pos)/10))+ 1/(10^(spn_pn_chan_conv_freq_SNR(mse_pos)/10));
      spn_ce_SNR(mse_pos) = 10*log10(1/temp1);
      mse_pos = mse_pos +1;
end
 spn_ce_SNR
dpn_ce_SNR

figure;
plot(abs(spn_h_pn_smooth_result));
title('PN信道估计结果');

figure;
subplot(1,2,1);
plot(abs(dpn_h_smooth_result));
title('双PN信道估计结果');
subplot(1,2,2);
plot(abs(spn_h_smooth_result));
title('单PN信道估计结果');

figure;
subplot(1,2,1)
semilogy(SNR,spn_mean_mse,'r*-');
title('单PN估计MSE');
subplot(1,2,2)
semilogy(SNR,dpn_mean_mse,'k*-');
title('双PN估计MSE');

figure;hold on;
plot(SNR,spn_pnrm_SNR,'r*-');
plot(SNR,dpn_pnrm_SNR,'k*-');
legend('单PN','双PN');
title('数据循环重构后的信噪比');

figure;hold on;
plot(SNR,spn_pn_chan_conv_freq_SNR,'r*-');
plot(SNR,dpn_pn_chan_conv_freq_SNR,'k*-');
legend('单PN','双PN');
title('信道估计误差的信噪比');