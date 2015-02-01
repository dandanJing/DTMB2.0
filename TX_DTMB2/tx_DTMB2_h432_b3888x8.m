%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
clear all,close all,clc

%%参数定义
PN_len = 255;  % PN 长度
PN_total_Len = 432; %帧头长度,前同步88，后同步89
PN_cyclic_Len = PN_total_Len - PN_len;%帧头中循环扩展的长度
PN_power = 3; %帧头幅度dB
FFT_Len = 3888*8; %帧体所需的FFT、IFFT长度
Frame_len = PN_total_Len + FFT_Len; %帧长
Super_Frame = 10; %超帧长度
Srrc_oversample = 1; %过采样率
Symbol_rate = 7.56e6; %符号速率
Sampling_rate = Symbol_rate * Srrc_oversample;%采样速率
QAM = 0;    %  0: 64QAM ,2:256APSK
BitPerSym = 6;
sim_num=100; %仿真的帧数

%%帧头信号
PN = PN_gen*1.975;

%%帧体
data_bef_map=zeros(1,FFT_Len*BitPerSym);
data_aft_map_tx1 = zeros(1,Frame_len*sim_num);

%%channel
max_delay = 10;
doppler_freq = 0;
isFading = 0;
channelFilter = get_multipath_channel(101,isFading,max_delay,doppler_freq,Sampling_rate,Symbol_rate);
channelFilter = 1;

%%data
start_pos = 1;
data_transfer = zeros(1,sim_num*FFT_Len);
data_start_pos = 1;
for i=1:sim_num
    data_x = randint(1,FFT_Len*BitPerSym);
    modtemp1=map64q(data_x); %%星座映射
    modtemp= modtemp1*3888*20;
    temp_t1=ifft(modtemp, FFT_Len);
    data_transfer(data_start_pos:data_start_pos+FFT_Len-1)=modtemp;
    data_start_pos = data_start_pos + FFT_Len;
    
    frm_len = Frame_len;
    data_aft_map_tx1(start_pos:start_pos+frm_len-1)=[PN temp_t1];
    start_pos = start_pos+frm_len;
end

Send_data_srrc_tx1 = data_aft_map_tx1;
Send_data_srrc_tx1_ch = filter(channelFilter,1,Send_data_srrc_tx1);%过信道

%%接收机
chan_len = 260;%信道长度
MAX_CHANNEL_LEN = PN_total_Len;
last_frame_tail = [0];
stateSrrcReceive = [];
h_prev1 = []; %前一帧信道估计
h_prev2 = [];  %前两帧信道估计
recover_data = zeros(1,sim_num*(Frame_len-PN_total_Len));
recover_data_pos = 1;
start_pos = 1;
h_off_thresh = 0.2; %根据前两帧信道估计当前帧时设置的阈值
 for i=1:sim_num
      Receive_data = Send_data_srrc_tx1_ch((i-1)*Frame_len*Srrc_oversample+1:i*Frame_len*Srrc_oversample);
      
      if(i < 3)  %前两帧进行信道估计，不迭代
          chan_in = Receive_data(1:PN_total_Len+chan_len);
          h_current = channel_estimate(chan_in, PN);
          h_pn_conv = channel_pn_conv(PN,h_current,chan_len);
          %保存状态
          h_prev2 = h_prev1;
          h_prev1 = h_current;
	      h_pn_conv_prv = h_pn_conv;
          continue;
      end
      
      close all;
      figure;
      plot(abs(Receive_data));
      title('当前帧数据');
      
      %%上一帧数据估计
      [sc_h1 sc_h2 sc_ha] = h_estimate_A(h_prev2, h_prev1, i, h_off_thresh);
      h_pn_conv = channel_pn_conv(PN,sc_h1,chan_len);
      last_frame_data =  Send_data_srrc_tx1_ch((i-2)*Frame_len+1:(i-1)*Frame_len+chan_len);
      last_frame_data(1:length(h_pn_conv_prv))= last_frame_data(1:length(h_pn_conv_prv))-h_pn_conv_prv;
      last_frame_data(Frame_len+(1:chan_len))= last_frame_data(Frame_len+(1:chan_len))-h_pn_conv(1:chan_len);
      
      iter_num = 2;
      for k = 1 : iter_num
          %%op1 根据前两帧进行信道粗时域冲激响应估计
          if(k==1)
				
			else
				[sc_h2 sc_ha] = h_estimate_B(h_prev2, h_prev1, i, h_off_thresh);
          end
          
          %%op2 减掉PN和拖尾
          h_pn_conv = channel_pn_conv(PN,sc_h2,chan_len);
          data_pn_rm = Receive_data;
          data_pn_rm = data_pn_rm(1:length(h_pn_conv)) - h_pn_conv;
          
          %%op3计算数据头和拖尾
         
      end
      
      %%op1 
      PN_R_time = Receive_data(1:2*PN_total_Len);
      fft_PN_R = fft(PN_R_time);
      fft_PN = fft(PN,2*PN_total_Len);
      H_F =  fft_PN_R./fft_PN;
      
      %平滑滤波
      freq_thres = 3.6;
      mmse_weight = 0.99;
      H_F(abs(fft_PN)<freq_thres)=0;
      h_coarse = ifft(H_F);
      h_mmse_filter = channel_mmse_filter(h_coarse, mmse_weight);
      
      id = 1;
      if(id == 1)
          figure;
          plot(abs(h_coarse));
          hold on;
          plot(abs(h_mmse_filter),'r');
          hold off;
          pause;
      end
      
      if(i == 1)
          continue;
      end
      %%op2 粗均衡，获得信号时域
      data_r = Receive_data(PN_total_Len+1:Frame_len);
      fft_data = fft(data_r);
      h = zeros(1,length(data_r));
      h(1:length(h_mmse_filter)) = h_mmse_filter;
      fft_channel = fft(h,length(data_r));
      data_freq = fft_data./(fft_channel);
      data_time_eq = ifft(data_freq);
      
      %%op3 PN重构,信道精估计
      [Y I] = max(abs(h_coarse));
      current_frame_data_head = conv(data_time_eq(1:PN_total_Len),h_coarse);
      current_frame_data_head = current_frame_data_head(I:end);
      current_frame_data_tail = conv(data_time_eq(length(data_r)-PN_total_Len+1:end),h_coarse);
      channel_len = length(h_coarse);
      PN_reshape = PN_R_time  ;
      %去除当前帧数据
      PN_reshape(PN_total_Len+1:2*PN_total_Len) = PN_reshape(PN_total_Len+1:2*PN_total_Len)-current_frame_data_head(1:PN_total_Len);
      %去除上一帧拖尾
      %PN_reshape(1:length(last_frame_tail)) = PN_reshape(1:length(last_frame_tail))-last_frame_tail;
    
      PN_receive_freq = fft(PN_reshape,2*PN_total_Len);
      PN_freq_datalen = fft(PN,2*PN_total_Len);
      channel_freq2 = PN_receive_freq./PN_freq_datalen;
      channel_estimate_time = ifft(channel_freq2);
      channel_estimate_freq = fft(channel_estimate_time,length(data_r));
      
      %%op4 数据循环重构
      data_reshape = data_r;
      data_reshape(1:PN_len)=data_reshape(1:PN_len)-PN_reshape(PN_len+1:2*PN_len)+current_frame_data_tail(PN_len+1:2*PN_len);
      data_reshape_freq = fft(data_reshape);
      data_estimate_freq = data_reshape_freq./channel_estimate_freq;
      data_estimate_time = ifft(data_estimate_freq);
      recover_data(start_pos:start_pos+length(data_estimate_time)-1) = data_estimate_time;
      start_pos = start_pos +  length(data_estimate_time);
      current_frame_data_tail = conv(data_estimate_time(length(data_estimate_time)-PN_total_Len+1:end),channel_estimate_time(1:2*PN_total_Len));
      last_frame_tail = current_frame_data_tail(PN_total_Len+1:2*PN_total_Len);
 end