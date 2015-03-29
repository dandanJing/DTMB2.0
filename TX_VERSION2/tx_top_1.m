%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48*8
clear all,close all,clc

debug = 0;
debug_multipath = 0;%定义是否考虑多径
debug_path_type = 16;%定义多径类型
SNR = [15:5:25];

%%参数定义
PN_len = 255;  % PN 长度
PN_total_len = 432; %帧头长度,前同步88，后同步89
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
BitPerSym = 6;
Srrc_oversample = 1; %过采样率
sim_num=1000; %仿真的帧数

%%帧体
data_bef_map=zeros(1,FFT_len*BitPerSym);

%%channel
max_delay = 10;
doppler_freq = 0;
isFading = 0;
channelFilter = 1;
if debug_multipath
    channelFilter = multipath_new(debug_path_type,1/7.56,1,0);
end

%%data
start_pos = 1;
data_transfer = zeros(1,sim_num*FFT_len);
data_start_pos = 1;

%%TPS value
load TPS_symbol.mat
for SNR_IN = SNR  %定义输入信噪比
    for i = 1:sim_num
        data_x = randi([0 1],1,FFT_len*BitPerSym);
        modtemp1=map64q(data_x); %%星座映射
        modtemp= modtemp1*3888*20;

        %%TPS
        tps_position =TPS_pos_gen(i,0);
        modtemp(tps_position)= tps_value;
        temp_t1=ifft(modtemp, FFT_len);
        data_transfer(data_start_pos:data_start_pos+FFT_len-1)=modtemp;
        data_start_pos = data_start_pos + FFT_len;

        PN_temp = PN_gen(i,0);
        frm_len = Frame_len;
        data_aft_map_tx1(start_pos:start_pos+frm_len-1)=[PN_temp temp_t1];
        start_pos = start_pos+frm_len;
    end
    Send_data_srrc_tx1 = data_aft_map_tx1;
    Send_data_srrc_tx1_ch_spn = filter(channelFilter,1,Send_data_srrc_tx1);%过信道
    Send_data_srrc_tx1_spn = awgn(Send_data_srrc_tx1_ch_spn,SNR_IN,'measured');

    matfilename = strcat('DTMB2_data_awgn_SNR',num2str(SNR_IN),'.mat');
    if debug_multipath
        matfilename = strcat('DTMB2_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
    end
    save(matfilename,'Send_data_srrc_tx1_spn','data_transfer');
end