%%channel estimation 
%%论文 1.Channel Estimation for the Chinese DTTB
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
%%方法2
clear all,close all,clc
debug = 0;
debug_tps = 1;
SNR = 20;

spn_mean_mse = zeros(1,length(SNR));
dpn_mean_mse = zeros(1,length(SNR));
mse_pos = 1;
for SNR_IN = SNR %定义输入信噪比
    
%%数据输入
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 1;%定义多径类型
if debug_multipath
    matfilename = strcat('DTMB_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
    load(matfilename);
else
    load strcat('DTMB_data_awgn_SNR',num2str(SNR_IN),'.mat');
end

debug_frame_32k_eq = 1;%定义数据帧均衡长度，1为32K， 0为FFT_Len

%%参数定义
PN_len = 255;  % PN 长度
PN_total_len = 432; %帧头长度,前同步88，后同步89
DPN_total_len = 1024;
DPN_len = 512;
load pn256_pn512.mat
PN_cyclic_Len = PN_total_len - PN_len;%帧头中循环扩展的长度
PN_power = 3; %帧头幅度dB
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
Super_Frame = 10; %超帧长度
Srrc_oversample = 1; %过采样率
Symbol_rate = 7.56e6; %符号速率
Sampling_rate = Symbol_rate * Srrc_oversample;%采样速率
QAM = 0;    %  0: 64QAM ,2:256APSK
BitPerSym = 6;
sim_num=1000; %仿真的帧数
iter_num = 2; %迭代次数

%%帧头信号
PN = PN_gen*1.975;

%%接收机
chan_len = 260;%信道长度
MAX_CHANNEL_LEN =PN_total_len;
last_frame_tail = [0];
stateSrrcReceive = [];
h_prev1 = []; %前一帧信道估计
h_prev2 = [];  %前两帧信道估计
recover_data = zeros(1,sim_num*(Frame_len-PN_total_len));
recover_data_pos = 1;
start_pos = 1;
h_off_thresh = 0.02; %根据前两帧信道估计当前帧时设置的阈值

debug_h_ave = 1;%判断是否对信道平滑
h_average = zeros(1,2048); %信道估计的平均结果
h_ave_alpha = 0.95;
h_start_ave_frame = 15;
h_start_iter_frame = 105;
h_average_thresh = 0.005;

channel_estimate_spn = zeros(sim_num,MAX_CHANNEL_LEN);
channel_estimate_dpn = zeros(sim_num,MAX_CHANNEL_LEN);
channel_estimate_temp =  zeros(sim_num,MAX_CHANNEL_LEN);
%%单PN信道估计
 for i=1:sim_num
      close all;
      Receive_data = Send_data_srrc_tx1_ch((i-1)*Frame_len+1:i*Frame_len);
      
      if(i < 3)  %前两帧进行信道估计，不迭代
          chan_in = Receive_data(1:PN_total_len+chan_len);
          h_current = channel_estimate(chan_in, PN, 2048, 0.1,debug);
          h_pn_conv = channel_pn_conv(PN,h_current,chan_len);
          %保存状态
          h_prev2 = h_prev1;
          h_prev1 = h_current;
	      h_pn_conv_prv = h_pn_conv;
          channel_estimate_spn(i,:) = h_current(1:PN_total_len);
          continue;
      end
       
      %%上一帧数据估计
      [sc_h1 sc_h2 sc_ha] = h_estimate_A(h_prev2, h_prev1, i, h_off_thresh);
      h_pn_conv = channel_pn_conv(PN,sc_h1,chan_len);
      last_frame_data =  Send_data_srrc_tx1_ch((i-2)*Frame_len+1:(i-1)*Frame_len);
      last_frame_data_tail = Send_data_srrc_tx1_ch((i-1)*Frame_len+(1:chan_len))- h_pn_conv(1:chan_len);
      last_frame_pn_tail = h_pn_conv_prv(PN_total_len+(1:chan_len));
      
      if debug_frame_32k_eq
           last_frame_data_tail_head = [ last_frame_data(PN_total_len+1:end),last_frame_data_tail ];
           last_frame_data_tail_head(1:chan_len) = last_frame_data_tail_head(1:chan_len)-last_frame_pn_tail;
           last_frame_ofdm_freq = fft(last_frame_data_tail_head, 32*1024);
           last_frame_h_freq = fft(sc_ha,32*1024);
      else
           last_frame_data_tail_head =  last_frame_data(PN_total_len+1:end);
           last_frame_data_tail_head(1:chan_len) = last_frame_data_tail_head(1:chan_len)-last_frame_pn_tail+ last_frame_data_tail;
           last_frame_ofdm_freq = fft(last_frame_data_tail_head);
           last_frame_h_freq = fft(sc_ha,FFT_len);
      end
       
     last_frame_ofdm_eq =  last_frame_ofdm_freq./last_frame_h_freq;
     last_frame_ofdm_eq_data = ifft(last_frame_ofdm_eq);
     last_frame_ofdm_eq_data =last_frame_ofdm_eq_data(1:FFT_len);
     fft_data = fft(last_frame_ofdm_eq_data);
     if debug_tps 
         channel_temp = channel_estimate_spn(i-1,:);
         channel_freq = fft(channel_temp, FFT_len);
         channel_freq(tps_position)= fft_data(tps_position)./tps_symbol;
         fft_data(tps_position)=tps_symbol;
     end
     recover_data((i-2)*FFT_len+1:(i-1)*FFT_len)=  fft_data;
     last_frame_ofdm_eq_data = ifft(fft_data);
      
      if debug || i== sim_num
          figure;
          plot(fft_data,'.b');
          title('均衡后的数据');
      end
      h_iter = sc_h1;
      last_frame_ofdm_eq_freq = fft(last_frame_ofdm_eq_data, 32*1024);
      last_frame_h_32k = fft(sc_ha,32*1024);
      last_frame_ofdm_h =  last_frame_ofdm_eq_freq.* last_frame_h_32k;
      last_frame_ofdm_h_conv = ifft(last_frame_ofdm_h);
      last_frame_data_tail = last_frame_ofdm_h_conv(FFT_len+1:FFT_len+chan_len);
      current_frame_pn = Receive_data(1:PN_total_len);
      current_frame_pn(1:chan_len)=current_frame_pn(1:chan_len)-last_frame_data_tail;
      for k = 1 : iter_num
          if k==1
            channel_estimate_thresh = 0.06;
          elseif k==2
            channel_estimate_thresh = 0.05;
          else
            channel_estimate_thresh = 0.05-0.01*(k-2);
            if 0.05-0.01*(k-2) < 0.01
                channel_estimate_thresh = 0.01;
            end
          end
          
          h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
          pn_recover = [current_frame_pn h_pn_conv(PN_total_len+(1:chan_len))];
          if k==1
              h_iter = channel_estimate(pn_recover,PN, 2048,channel_estimate_thresh,debug);
          else
              h_iter = channel_estimate_A(pn_recover,PN, 2048,debug);
          end
      end
      
      h_iter_old = sc_h1;
      
      h_iter_old = h_iter;
      h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
      pn_recover = [current_frame_pn h_pn_conv(PN_total_len+(1:chan_len))];
      h_iter = channel_estimate_A(pn_recover,PN, 2048,debug);
      
      if debug_h_ave && i > h_start_ave_frame 
            h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
            pn_recover = [current_frame_pn h_pn_conv(PN_total_len+(1:chan_len))];
            h_temp = channel_estimate_B(pn_recover,PN, 2048,debug);
            channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
            if i == h_start_iter_frame
                h_average(1:PN_total_len) = mean(channel_estimate_temp(h_start_ave_frame+1:h_start_iter_frame,:));
                h_average = channel_denoise(h_average,h_average_thresh);
                h_iter = h_average;
            elseif i > h_start_iter_frame
                h_temp_1 = h_temp;
                h_temp_1(PN_total_len+1:end)=0;
                h_iter = h_ave_alpha*h_average+(1-h_ave_alpha)*h_temp_1;
                h_iter = channel_denoise(h_iter,h_average_thresh);
                h_average = h_iter;
            end
      end
      
      if debug
      figure;
      subplot(1,2,1);
      plot(abs(h_iter_old(1:PN_total_len)),'r');
      title('初始信道估计');
       subplot(1,2,2);
       plot(abs(h_iter(1:PN_total_len)),'b');
       title('迭代信道估计');
       if debug
             pause;
        end
      end
      
      chan_len = min(chan_len_estimate(h_iter),MAX_CHANNEL_LEN);
      h_prev2 = h_prev1;
      h_prev1 = h_iter;
	  h_pn_conv_prv = channel_pn_conv(PN,h_iter,chan_len);
      channel_estimate_spn(i,:) = h_iter(1:PN_total_len);
 end
 
 %%真实多径信道与单PN估计结果比较
 if debug_multipath
     max_delay = 10;
    doppler_freq = 0;
    isFading = 0;
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
    channel_real = zeros(1,PN_total_len);
    channel_real(1:length(channelFilter)) = channelFilter;
    figure;
    subplot(1,2,1);
    plot(abs(channel_real));
    title('真实多径信道');
    subplot(1,2,2);
    plot(abs( h_iter(1:PN_total_len)),'b');
    title('迭代信道估计');   
 else
     channel_real = zeros(1,PN_total_len);
     channel_real(1) = 1;
 end
 
  %% 双PN信道估计
 frame_test = DPN_total_len + FFT_len;
 for i=1:sim_num-1
      Receive_data = Send_data_srrc_tx1_ch2((i-1)*frame_test+(1:frame_test));
      pn_test = Receive_data(DPN_len+(1:DPN_len));
      coeff = 6.4779e+04;
      pn_test = pn_test ./ coeff;
      pn512_fft = fft(pn_test);
      dpn_h_freq =  pn512_fft./ pn512;
      dpn_h_time = ifft(dpn_h_freq);
      chan_len_test = chan_len_estimate(dpn_h_time);
      channel_estimate_dpn(i,1:chan_len_test)=dpn_h_time(1:chan_len_test);
 end
 
 %%估计误差
 num = sim_num-9;
 mean_pos = h_start_ave_frame+1:num;
 dpn_channel_mean = mean(channel_estimate_dpn(mean_pos,:));
 dpn_mean_off = dpn_channel_mean -  channel_real;
 dpn_mean_mse(mse_pos) = norm(dpn_mean_off)/norm(channel_real);
 
 spn_channel_mean = mean(channel_estimate_spn(mean_pos,:));
 spn_mean_off = spn_channel_mean -  channel_real;
 spn_mean_mse(mse_pos) = norm(spn_mean_off)/norm(channel_real);

 temp_mean = mean(channel_estimate_temp(mean_pos,:));
 temp_off = temp_mean -  channel_real;
 temp_mse(mse_pos) =  norm(temp_off)/norm(channel_real);
 
 if debug || i== sim_num-1 
        figure;
        subplot(1,2,1);
        plot(abs(channel_real));
        title('真实多径信道');
        subplot(1,2,2);
        plot(abs(dpn_channel_mean));
        title('双PN信道平均估计结果');
        figure;
        subplot(1,2,1);
        plot(abs(channel_real));
        title('真实多径信道');
        subplot(1,2,2);
        plot(abs(spn_channel_mean));
        title('单PN信道平均估计结果');
        figure;
        subplot(1,2,1);
        plot(abs(channel_real));
        title('真实多径信道');
        subplot(1,2,2);
        plot(abs(temp_mean));
        title('temp平均估计结果');
        if debug
            pause;
        end
 end 
 mse_pos = mse_pos + 1;
end

figure;
subplot(1,3,1)
semilogy(SNR(2:end),spn_mean_mse(2:end),'r*-');
title('单PN估计MSE');
subplot(1,3,2)
semilogy(SNR,dpn_mean_mse,'k^-');
title('双PN信道平均估计MSE');
subplot(1,3,3)
semilogy(SNR,temp_mse,'k^-');
title('temp信道平均估计MSE');