%%channel estimation 
%%论文 
%%数据干扰消除
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
%%方法2
clear all,close all,clc
debug = 0;
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

    define_eq_total = 1;
    %%单PN信道估计
     for i=1:1
          close all;
          Receive_data = Send_data_srrc_tx1_ch((i-1)*Frame_len+1:i*Frame_len);

          chan_in = Receive_data(1:2*PN_total_len);
          h_current = channel_estimate(chan_in, PN, 2048, 0.1,debug);
          for k = 1:iter_num
              if define_eq_total
              else
              end
          end
     end
end
 figure;
 plot(abs(h_current));
 