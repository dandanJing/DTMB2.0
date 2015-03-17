%%channel estimation 
%%论文 
%%数据干扰消除
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
%%方法2
clear all,close all,clc
debug = 0;
debug_iter =0;
debug_tps = 1;
debug_eq_total = 0;%1为 frame_len  0为fft_len
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 8;%定义多径类型
SNR = [20];

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
iter_num = 3; %迭代次数
MAX_CHANNEL_LEN = PN_total_len;

%%帧头信号
PN = PN_gen*1.975;
temp = ifft(pn512);
DPN = temp*sqrt(var(PN)/var(temp));

%%信道
if debug_multipath
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
    channel_real = zeros(1,PN_total_len);
    channel_real(1:length(channelFilter)) = channelFilter;
else
    channel_real = zeros(1,PN_total_len);
    channel_real(1) = 1;
end

spn_mean_mse = zeros(1,length(SNR));
dpn_mean_mse = zeros(1,length(SNR));
spn_pnrm_SNR = zeros(1,length(SNR));
dpn_pnrm_SNR = zeros(1,length(SNR));
spn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));
dpn_pn_chan_conv_freq_SNR = zeros(1,length(SNR));

mse_pos = 1;
for SNR_IN = SNR %定义输入信噪比
   
    %%数据输入
    if debug_multipath
        matfilename = strcat('DTMB_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
        load(matfilename);
    else
        load strcat('DTMB_data_awgn_SNR',num2str(SNR_IN),'.mat');
    end
    
    %%接收机
    chan_len = PN_total_len;%信道长度
    last_frame_tail = [0];
    stateSrrcReceive = [];
    h_prev1 = []; %前一帧信道估计
    h_prev2 = [];  %前两帧信道估计
    recover_data = zeros(1,sim_num*(Frame_len-PN_total_len));
    recover_data_pos = 1;
    start_pos = 1;
    h_off_thresh = 0.02; %根据前两帧信道估计当前帧时设置的阈值

    h_start_frame_num = 105;
    h_ave_frame_num = 50;
    h_ave_frame_num_dpn = 100;
    h_start_ave_frame_num = 5;
    h_average_thresh = 0.005;
    h_denoise_thresh = 0.008;
    
    channel_estimate_spn = zeros(sim_num,MAX_CHANNEL_LEN);
    channel_estimate_dpn = zeros(sim_num,MAX_CHANNEL_LEN);
    channel_estimate_temp =  zeros(sim_num,MAX_CHANNEL_LEN);

    last_frame_data_conv = zeros(1,2048);
    start_pos = 0;
    %%单PN信道估计
     for i=1:sim_num-1
          close all;
          figure;
          plot(abs(channel_real));
          title('真实多径信道');

          last_frame_start_pos = start_pos - Frame_len;
          if mod(i,Super_Frame)==1
              start_pos = start_pos + PN_total_len;
          end
          
          %当前帧数据
          Receive_data = Send_data_srrc_tx1_ch(start_pos+(1:Frame_len));

          %迭代信道估计
          if mod(i,Super_Frame)==1 %超帧第一帧时不迭代
              chan_in = Receive_data(1:PN_total_len);
              pn_freq = fft(chan_in);
              h_freq = pn_freq./fft(PN);
              h_current = ifft(h_freq);
              h_old = h_current;
              chan_len = min(chan_len_estimate(h_current),MAX_CHANNEL_LEN);
              if chan_len > 300
                  chan_len = min(chan_len_estimate_1(h_current),MAX_CHANNEL_LEN);
              end
              h_current(chan_len+1:end)=0;
              h_current = channel_denoise(h_current,h_denoise_thresh);
              if debug
                chan_off = h_current(1:PN_total_len)-channel_real;
                mse = norm(chan_off)/norm(channel_real)
                chan_len
                figure;
                plot(abs(h_old(1:PN_total_len)));
                title('超帧估计初始结果');
                figure;
                plot(abs(h_current(1:PN_total_len)));
                title('超帧头估计结果');
                pause;
              end
              channel_estimate_temp(i,:) = h_current(1:PN_total_len);
              if i >= h_ave_frame_num+h_start_ave_frame_num
                  h_iter = mean(channel_estimate_temp(i-h_ave_frame_num+1:i,:));
              else
                  h_iter = h_current(1:PN_total_len);
              end
          else
%               if i > 2
%                   [sc_h1 sc_h2] = h_estimate_A(h_prev2, h_prev1, i, h_off_thresh);
%                   sc_ha = (sc_h1 + sc_h2)/2;
%               else
%                   sc_ha = h_prev1;
%                   sc_h1 = h_prev1;
%               end
%               h_iter = sc_ha;
              for k = 1:iter_num
                  if debug_eq_total
                     eq_in_data = Send_data_srrc_tx1_ch(start_pos+PN_total_len+(1:Frame_len));
                     eq_data_freq = fft(eq_in_data);
                     h_freq = fft( h_iter,Frame_len);
                     eq_data_freq = (eq_data_freq./h_freq);
                     eq_data_freq(abs(h_freq)<h_off_thresh) =  0;
                     eq_data_time = ifft(eq_data_freq);
                     h_freq_fftlen =  fft( h_iter,2*PN_total_len);
                     fft_data_freq = fft(eq_data_time(1:PN_total_len),2*PN_total_len);
                     data_chan_conv = ifft(fft_data_freq.* h_freq_fftlen);
                     chan_in = Receive_data(1:PN_total_len+chan_len);
                     chan_in(1:chan_len)= chan_in(1:chan_len)-last_frame_data_conv(PN_total_len+(1:chan_len));
                     chan_in(PN_total_len+(1:chan_len))=chan_in(PN_total_len+(1:chan_len))-data_chan_conv(1:chan_len);
                     h_temp = channel_estimate_B(chan_in,PN, 2048,debug);
                     h_temp(chan_len+1:end)=0;
                     channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
                     if i < h_ave_frame_num+h_start_ave_frame_num
                         h_iter = channel_estimate_temp(1,:);
                         h_temp = channel_denoise(h_temp,h_denoise_thresh);
                         channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
                     else
                         h_iter = mean(channel_estimate_temp(i-h_ave_frame_num+1:i,:));
                         if chan_len > 250
                              h_iter = mean(channel_estimate_temp(h_start_ave_frame_num+1:i,:));
                         end
                     end
                  else
                      h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
                      eq_in_data = Receive_data(PN_total_len+1:end);
                      data_tail =  Send_data_srrc_tx1_ch(start_pos+Frame_len+(1:chan_len))-h_pn_conv(1:chan_len);
                      eq_in_data(1:chan_len)=eq_in_data(1:chan_len)-h_pn_conv(PN_total_len+(1:chan_len))+data_tail;
                      eq_in_data_freq = fft(eq_in_data);
                      h_freq = fft( h_iter,FFT_len);
                      eq_data_freq = eq_in_data_freq./h_freq;
                      eq_data_freq(abs(h_freq)<h_off_thresh) = 0;
                      eq_data_time = ifft(eq_data_freq);
                      h_freq_fftlen =  fft( h_iter,2*PN_total_len);
                      fft_data_freq = fft(eq_data_time(1:PN_total_len),2*PN_total_len);
                      data_chan_conv = ifft(fft_data_freq.* h_freq_fftlen);
                      chan_in = Receive_data(1:PN_total_len+chan_len);
                      chan_in(1:chan_len)= chan_in(1:chan_len)-last_frame_data_conv(PN_total_len+(1:chan_len));
                      chan_in(PN_total_len+(1:chan_len))=chan_in(PN_total_len+(1:chan_len))-data_chan_conv(1:chan_len);
                      h_temp = channel_estimate_B(chan_in,PN, 2048,debug);
                      h_temp(chan_len+1:end)=0;
                     channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
                     if i < h_ave_frame_num+h_start_ave_frame_num
                         h_iter = channel_estimate_temp(1,:);
                         h_temp = channel_denoise(h_temp,h_denoise_thresh);
                         channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
                     else
                         h_iter = mean(channel_estimate_temp(i-h_ave_frame_num+1:i,:));
                         if chan_len > 250
                              h_iter = mean(channel_estimate_temp(h_start_ave_frame_num+1:i,:));
                         end
                     end
                  end
              end
          end
          chan_len = min(chan_len_estimate(h_iter),MAX_CHANNEL_LEN);
          h_iter(chan_len+1:end)=0;
          if chan_len > 250
              h_iter = channel_denoise(h_iter,h_denoise_thresh);
          end
          channel_estimate_spn(i,:) = h_iter(1:PN_total_len);
          h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
          
          h_sa = h_iter(1:PN_total_len);
          h_freq = fft(h_sa,FFT_len);
          fft_data = Receive_data(PN_total_len+1:end);
          fft_data_tail =  Send_data_srrc_tx1_ch(start_pos+Frame_len+(1:chan_len))-h_pn_conv(1:chan_len);
          fft_PN_tail = h_pn_conv(PN_total_len+(1:chan_len));
          fft_data(1:chan_len)=fft_data(1:chan_len)-fft_PN_tail+fft_data_tail;
          fft_data_freq = fft(fft_data);
          
          ofdm_eq =  fft_data_freq./h_freq;
          max_h = max(abs(h_freq));
          ofdm_eq(abs(h_freq)<max_h*h_average_thresh)=0;
          frame_eq_data_time = ifft(ofdm_eq);
      
          data_time_tail = frame_eq_data_time(FFT_len-PN_total_len+1:FFT_len);
          data_h_freq_temp = fft(h_sa, 2*PN_total_len);
          data_freq_temp = fft(data_time_tail, 2*PN_total_len);
          last_frame_data_conv = ifft(data_freq_temp.*data_h_freq_temp);
          
          %数据重构后的时域数据恢复结果
          spn_rcov_channel_data_time((i-1)*FFT_len+(1:FFT_len))= fft_data;
          %数据重构后的频域数据恢复结果
          spn_ch_data_conv_freq((i-1)*FFT_len+(1:FFT_len))=h_freq.*data_transfer((i-1)*FFT_len+(1:FFT_len));
          
          if debug||i==sim_num-1||(debug_iter&&i > h_ave_frame_num+h_start_ave_frame_num-2)
              i
              chan_off = h_iter(1:PN_total_len)-channel_real;
              mse = norm(chan_off)/norm(channel_real)
              figure;
              plot(ofdm_eq,'.k');
              title('单PN均衡后的数据');
              figure;
              plot(abs(h_iter));
              title('迭代估计结果');
              if ~(i==sim_num-1)
                pause
              end
          end
          h_prev2 = h_prev1;
          h_prev1 = h_iter;
          start_pos = start_pos + Frame_len;
     end
     
      %% 双PN信道估计
      %真实数据与真实信道的卷积结果
     rcov_channel_real_data_time = zeros(1,FFT_len*(sim_num-1));
     rcov_channel_real_data_freq = zeros(1,FFT_len*(sim_num-1));

     %双PN去除PN拖尾后的数据恢复结果
     dpn_rcov_channel_data_time = zeros(1,FFT_len*(sim_num-1));
     %双PN信道估计结果与真实传输数据的卷积结果
     dpn_ch_data_conv_freq = zeros(1,FFT_len*(sim_num-1));

     channel_real_freq = fft(channel_real,FFT_len);
     frame_test = DPN_total_len + FFT_len;
     coeff = 6.4779e+04;
     for i=1:sim_num-1
          Receive_data = Send_data_srrc_tx1_ch2((i-1)*frame_test+(1:frame_test));
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
          chan_len_dpn = min(chan_len_estimate(dpn_h_time)+20,MAX_CHANNEL_LEN);
          dpn_h_time(chan_len_dpn+1:end)=0;
          channel_estimate_dpn(i,1:PN_total_len)=dpn_h_time(1:PN_total_len);

          temp_PN = DPN;
          temp_PN_len = DPN_len;
          temp_total_len = DPN_total_len;

          frame_test =temp_total_len + FFT_len;
          Receive_data = Send_data_srrc_tx1_ch2((i-1)*frame_test+(1:frame_test));

          data_real =  data_transfer((i-1)*FFT_len+(1:FFT_len)).*channel_real_freq;
          rcov_channel_real_data_freq((i-1)*FFT_len+(1:FFT_len))= data_real;
          rcov_channel_real_data_time((i-1)*FFT_len+(1:FFT_len))= ifft(data_real);

          %%数据恢复
          if i >=  h_ave_frame_num_dpn
              dpn_chan_mean = mean(channel_estimate_dpn(i-  h_ave_frame_num_dpn+1:i,:));
          elseif i==1
              dpn_chan_mean = channel_estimate_dpn(1,:);
          else
              dpn_chan_mean = mean(channel_estimate_dpn(1:i,:));
          end
          chan_len_dpn = min(chan_len_estimate(dpn_chan_mean),MAX_CHANNEL_LEN)+10;
          dpn_chan_mean(chan_len_dpn+1:end)=0;
          dpn_h_freq_frame = fft(dpn_chan_mean,FFT_len);

          pn_conv = channel_pn_conv( temp_PN, dpn_chan_mean,chan_len_dpn);
          frame_tail =  Send_data_srrc_tx1_ch2(i*frame_test+(1:chan_len_dpn))-pn_conv(1:chan_len_dpn);
          frame_data =  Receive_data(temp_total_len+1:end);
          frame_data(1:chan_len_dpn) = frame_data(1:chan_len_dpn)-pn_conv(temp_PN_len+(1:chan_len_dpn))+frame_tail;

          dpn_rcov_channel_data_time((i-1)*FFT_len+(1:FFT_len))=  frame_data;
          dpn_ch_data_conv_freq((i-1)*FFT_len+(1:FFT_len))=data_transfer((i-1)*FFT_len+(1:FFT_len)).*dpn_h_freq_frame;

          fft_data_dpn = fft(frame_data);
          fft_data_dpn = fft_data_dpn./dpn_h_freq_frame;
          recover_data_dpn((i-1)*FFT_len+1:(i)*FFT_len)=  fft_data_dpn;
     end
      figure;
      plot(fft_data_dpn,'.k');
      title('双PN数据均衡结果');
      
      %%估计误差
     num = sim_num-9;
     mean_pos =  h_start_frame_num+1:num;
     dpn_channel_mean = mean(channel_estimate_dpn(mean_pos,:));
     chan_len = min(chan_len_estimate(dpn_channel_mean)+10,MAX_CHANNEL_LEN);
     dpn_channel_mean(chan_len+1:end)=0;
     dpn_mean_off = dpn_channel_mean -  channel_real;
     dpn_mean_mse(mse_pos) = norm(dpn_mean_off)/norm(channel_real);

     spn_channel_mean = mean(channel_estimate_spn(mean_pos,:));
      spn_channel_mean(chan_len+1:end)=0;
     spn_mean_off = spn_channel_mean -  channel_real;
     spn_mean_mse(mse_pos) = norm(spn_mean_off)/norm(channel_real);

     temp_mean = mean(channel_estimate_temp(mean_pos,:));
     chan_len = min(chan_len_estimate(temp_mean),MAX_CHANNEL_LEN);
     temp_mean(chan_len+1:end)=0;
     temp_off = temp_mean -  channel_real;
     temp_mse(mse_pos) =  norm(temp_off)/norm(channel_real);
 
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
     start_pos = FFT_len*h_start_frame_num+1;
     end_pos = FFT_len*(sim_num-9);
     dpn_pnrm_SNR(mse_pos) = estimate_SNR(dpn_rcov_channel_data_time(start_pos:end_pos),rcov_channel_real_data_time(start_pos:end_pos))
     spn_pnrm_SNR(mse_pos) = estimate_SNR(spn_rcov_channel_data_time(start_pos:end_pos),rcov_channel_real_data_time(start_pos:end_pos))
     spn_pn_chan_conv_freq_SNR(mse_pos) = estimate_SNR(spn_ch_data_conv_freq(start_pos:end_pos),rcov_channel_real_data_freq(start_pos:end_pos))
     dpn_pn_chan_conv_freq_SNR(mse_pos) = estimate_SNR(dpn_ch_data_conv_freq(start_pos:end_pos),rcov_channel_real_data_freq(start_pos:end_pos))
     mse_pos = mse_pos + 1;
end
figure;
subplot(1,2,1)
semilogy(SNR,temp_mse,'r*-');
title('单PN估计MSE');
subplot(1,2,2)
semilogy(SNR,dpn_mean_mse,'k*-');
title('双PN估计MSE');

figure;hold on;
plot(SNR,spn_pnrm_SNR,'r*-');
plot(SNR,dpn_pnrm_SNR,'k*-');
legend('单PN','双PN');
title('数据循环重构后的信噪比');
hold off;

figure;hold on;
plot(SNR,spn_pn_chan_conv_freq_SNR,'r*-');
plot(SNR,dpn_pn_chan_conv_freq_SNR,'k*-');
legend('单PN','双PN');
title('信道估计误差的信噪比'); 
 