%%channel estimation 
%%论文 
%%数据干扰消除
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
%%方法2
clear all,close all,clc
debug = 1;
debug_tps = 1;
SNR = [25];

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
    iter_num = 3; %迭代次数

    %%帧头信号
    PN = PN_gen*1.975;

    %%接收机
    chan_len = PN_total_len;%信道长度
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

    define_eq_total = 0;
    last_frame_data_conv = zeros(1,2048);
    %%单PN信道估计
     for i=1:sim_num
          close all;
          Receive_data = Send_data_srrc_tx1_ch((i-1)*Frame_len+1:i*Frame_len);

          chan_in = Receive_data(1:PN_total_len+chan_len);
          chan_in(1:chan_len)= chan_in(1:chan_len)-last_frame_data_conv(PN_total_len+(1:chan_len));
          h_current = channel_estimate(chan_in, PN, 2048, 0.1,0);
          chan_len = min(chan_len_estimate(h_current),MAX_CHANNEL_LEN);
          h_current(chan_len+1:end)=0;
          h_pn_conv = channel_pn_conv(PN,h_current,chan_len);
          h_iter = h_current;
          for k = 1:iter_num
              if define_eq_total
                 eq_in_data = Send_data_srrc_tx1_ch((i-1)*Frame_len+PN_total_len+(1:Frame_len));
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
                 if i < h_start_iter_frame
                     if k == iter_num
                        h_iter = channel_estimate_A(chan_in, PN, 2048, 0);
                     else
                        h_iter = channel_estimate(chan_in, PN, 2048, 0.1, 0);
                     end
                 else
                     h_temp = channel_estimate_B(chan_in,PN, 2048,debug);
                     chan_len = min(chan_len_estimate(h_temp),MAX_CHANNEL_LEN);
                     h_temp(chan_len+1:end)=0;
                     channel_estimate_temp(i,:) = h_temp(1:PN_total_len);
                 end
                 
                 chan_len = min(chan_len_estimate(h_iter),MAX_CHANNEL_LEN);
                 h_iter(chan_len+1:end)=0;
                 h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
              else
                  eq_in_data = Receive_data(PN_total_len+1:end);
                  data_tail = Send_data_srrc_tx1_ch((i)*Frame_len+(1:chan_len))-h_pn_conv(1:chan_len);
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
                if k == iter_num
                    h_iter = channel_estimate_A(chan_in, PN, 2048, 0);
                 else
                    h_iter = channel_estimate(chan_in, PN, 2048, 0.1, 0);
                 end
                 chan_len = min(chan_len_estimate(h_iter),MAX_CHANNEL_LEN);
                 h_iter(chan_len+1:end)=0;
                 h_pn_conv = channel_pn_conv(PN,h_iter,chan_len);
              end
          end
          channel_estimate_spn(i,:) = h_iter(1:PN_total_len);
          if debug||i==sim_num
              figure;
              plot(abs(h_iter));
              title('迭代估计结果');
              if debug
                pause
              end
          end
     end
end
 
 