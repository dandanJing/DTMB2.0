%%iterative channel estimation   
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
clear all,close all,clc

debug = 0;
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 16;%定义多径类型
SNR = [20];
for SNR_IN = SNR  %定义输入信噪比

%%参数定义
PN_len = 255;  % PN 长度
PN_total_len = 432; %帧头长度,前同步88，后同步89
DPN_total_len = 1024;
PN_cyclic_Len = PN_total_len - PN_len;%帧头中循环扩展的长度
PN_power = 3; %帧头幅度dB
FFT_len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_len + FFT_len; %帧长
Super_Frame = 100; %超帧长度
Srrc_oversample = 1; %过采样率
Symbol_rate = 7.56e6; %符号速率
Sampling_rate = Symbol_rate * Srrc_oversample;%采样速率
QAM = 0;    %  0: 64QAM ,2:256APSK
BitPerSym = 6;
sim_num=100; %仿真的帧数

load pn256_pn512.mat
%%帧头信号
temp = ifft(pn512);
DPN = temp*6.4779*10^4;

%%帧体
data_bef_map=zeros(1,FFT_len*BitPerSym);
data_aft_map_tx1 = zeros(1,Frame_len*sim_num);

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
start_pos_1 = 1;
data_transfer = zeros(1,sim_num*FFT_len);
data_start_pos = 1;
for i=1:sim_num
    data_x = randi([0 1],1,FFT_len*BitPerSym);
    modtemp1=map64q(data_x); %%星座映射
    modtemp= modtemp1*3888*20;
    
    %%TPS
    [tps_position tps_symbol]=TPS_gen(i,0);
    modtemp(tps_position)= tps_symbol*1.975;
    temp_t1=ifft(modtemp, FFT_len);
    data_transfer(data_start_pos:data_start_pos+FFT_len-1)=modtemp;
    data_start_pos = data_start_pos + FFT_len;
    
    PN_temp = PN_gen(i,0)*1.975;
    frm_len = Frame_len;
    data_aft_map_tx1(start_pos:start_pos+frm_len-1)=[PN_temp temp_t1];
    start_pos = start_pos+frm_len;
    if debug
        close all;
        i
        figure;
        pn_freq = fft(PN_temp);
        min_value =  min(abs(pn_freq))
        plot(abs(pn_freq));
        title('当前帧PN频域响应');
    end
    
    frm_len_1 = FFT_len+DPN_total_len;
    data_aft_map_tx2(start_pos_1:start_pos_1+frm_len_1-1)=[DPN DPN temp_t1];
    start_pos_1 = start_pos_1+frm_len_1;
end

Send_data_srrc_tx1 = data_aft_map_tx1;
Send_data_srrc_tx1_ch_spn = filter(channelFilter,1,Send_data_srrc_tx1);%过信道
Send_data_srrc_tx1_spn = awgn(Send_data_srrc_tx1_ch_spn,SNR_IN,'measured');

Send_data_srrc_tx2 = data_aft_map_tx2;
Send_data_srrc_tx1_ch_dpn = filter(channelFilter,1,Send_data_srrc_tx2);%过信道
Send_data_srrc_tx1_dpn = awgn(Send_data_srrc_tx1_ch_dpn,SNR_IN,'measured');

matfilename = strcat('DTMB2_data_awgn_SNR',num2str(SNR_IN),'.mat');
if debug_multipath
    matfilename = strcat('DTMB2_data_multipath_new',num2str(debug_path_type),'SNR',num2str(SNR_IN),'.mat');
end
save(matfilename,'Send_data_srrc_tx1_spn','Send_data_srrc_tx1_ch_spn','Send_data_srrc_tx1_ch_dpn','Send_data_srrc_tx1_dpn','data_transfer',...
'Super_Frame');

end