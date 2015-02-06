%%iterative channel estimation    
%%DTMB2.0数据发送 帧头432，帧体3888*8，TPS 48,64QAM
clear all,close all,clc

debug = 1;
debug_multipath = 1;%定义是否考虑多径
debug_path_type = 102;%定义多径类型
SNR_IN = 20;%定义输入信噪比

%%参数定义
PN_len = 255;  % PN 长度
PN_total_len = 432; %帧头长度,前同步88，后同步89
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

%%帧头信号
PN = PN_gen*1.975;

%%帧体
data_bef_map=zeros(1,FFT_len*BitPerSym);
data_aft_map_tx1 = zeros(1,Frame_len*sim_num);

%%channel
max_delay = 10;
doppler_freq = 0;
isFading = 0;
channelFilter = 1;
if debug_multipath
    channelFilter = get_multipath_channel(debug_path_type,isFading,max_delay,doppler_freq,Sampling_rate,Symbol_rate);
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
    temp_t1=ifft(modtemp, FFT_len);
    data_transfer(data_start_pos:data_start_pos+FFT_len-1)=modtemp;
    data_start_pos = data_start_pos + FFT_len;
    
    frm_len = Frame_len;
    data_aft_map_tx1(start_pos:start_pos+frm_len-1)=[PN temp_t1];
    start_pos = start_pos+frm_len;
    
    frm_len_1 = Frame_len+PN_total_len;
    data_aft_map_tx2(start_pos_1:start_pos_1+frm_len_1-1)=[PN PN temp_t1];
    start_pos_1 = start_pos_1+frm_len_1;
end

Send_data_srrc_tx1 = data_aft_map_tx1;
Send_data_srrc_tx1_ch1 = filter(channelFilter,1,Send_data_srrc_tx1);%过信道
Send_data_srrc_tx1_ch = awgn(Send_data_srrc_tx1_ch1,SNR_IN,'measured');
SNR_Old = estimate_SNR(Send_data_srrc_tx1_ch,Send_data_srrc_tx1_ch1)

Send_data_srrc_tx2 = data_aft_map_tx2;
Send_data_srrc_tx1_ch2 = filter(channelFilter,1,Send_data_srrc_tx2);%过信道
Send_data_srrc_tx1_ch2 = awgn(Send_data_srrc_tx1_ch2,SNR_IN,'measured');

%%接收机
chan_len = 260;%信道长度
MAX_CHANNEL_LEN =PN_total_len;
last_frame_tail = [0];
stateSrrcReceive = [];
recover_data = zeros(1,sim_num*(Frame_len-PN_total_len));
recover_data_pos = 1;
start_pos = 1;
h_off_thresh = 0.1; %根据前两帧信道估计当前帧时设置的阈值
y_n_pre = zeros(1,chan_len);
channel_estimate_1 = zeros(sim_num,MAX_CHANNEL_LEN);
channel_estimate_2 = zeros(sim_num,MAX_CHANNEL_LEN);
 for i=1:sim_num-1
      Receive_data = Send_data_srrc_tx1_ch((i-1)*Frame_len+PN_total_len+(1:Frame_len));
      
      if(i < 3)  %前两帧进行信道估计，不迭代
          continue;
      end
      
      close all;
      pn_prefix =  Send_data_srrc_tx1_ch((i-1)*Frame_len+(1:PN_total_len));
      pn_prefix(1:chan_len) = pn_prefix(1:chan_len)-y_n_pre;
      z_n = zeros(1,2*PN_total_len);
      z_n(1:PN_total_len+chan_len)=[pn_prefix,Receive_data(1:chan_len)];
      iter_num = 2;
      h_iter = 0;
      alpha = 0.2;
      for k = 1:iter_num
          %%op1 channel estimate
          h_temp = channel_estimate(z_n, PN, 2*PN_total_len);
          h_temp = h_temp(1:PN_total_len);
          if k==1
              h_iter = h_temp;
          else
              h_iter = alpha*h_iter+(1-alpha)*h_temp;
          end
          
          if debug || (i== sim_num-1 && k == iter_num)
            figure;
            subplot(1,2,1);
            plot(abs(h_iter),'r');
            title('迭代估计结果');
            subplot(1,2,2);
            plot(abs(h_temp),'b');
           title('当前估计结果');
            hold off;
          end
          
          %%op2 data equalization
          frame_freq = fft( Receive_data,Frame_len);
          h_freq = fft(h_iter,Frame_len);
          frame_freq_eq = frame_freq./h_freq;
          frame_freq_eq(abs(h_freq)<h_off_thresh) =  0;
          frame_eq = ifft(frame_freq_eq);
          if debug || (i== sim_num-1 && k == iter_num)
              figure;
              plot(abs(h_freq));
              title('信道频域响应');
          end
          
        %% op3 data and channel convolution
         frame_ofdm_data = frame_eq(1:FFT_len);
         frame_ofdm_freq = fft(frame_ofdm_data, Frame_len);
         h_freq = fft(h_iter, Frame_len);
         y_n = ifft( frame_ofdm_freq.*h_freq);
         frame_recover_data =  fft(frame_ofdm_data, FFT_len);
         if debug || (i== sim_num-1 && k == iter_num)
            figure;
            plot(frame_recover_data,'.');
            title('数据均衡结果');
            if debug
                pause;
            end
         end
         
         z_n = zeros(1,2*PN_total_len);
         pn_tail_rm_data = Receive_data(1:chan_len)- y_n(1:chan_len);
         z_n(1:PN_total_len+chan_len)=[pn_prefix, pn_tail_rm_data ];
      end
       channel_estimate_1(i,1:length(h_iter))=h_iter;
       recover_data((i-1)*FFT_len+1:(i)*FFT_len)=  frame_recover_data;
       chan_len = min(chan_len_estimate(h_iter),MAX_CHANNEL_LEN);
       y_n_pre = y_n(FFT_len+(1:chan_len));
      
 end
 
 %% dual pn estimate 
 frame_test = PN_total_len + Frame_len;
 for i=1:sim_num-1
      Receive_data = Send_data_srrc_tx1_ch2((i-1)*frame_test+(1:frame_test));
      pn_test = Receive_data(PN_total_len+(1:PN_total_len));
      h_estimate = channel_estimate(pn_test,PN,PN_total_len);
      chan_len_test = chan_len_estimate( h_estimate);
      channel_estimate_2(i,1:chan_len_test)=h_estimate(1:chan_len_test);
 end
 
  figure;
 subplot(1,2,1);
 plot(abs(channel_estimate_1(sim_num-5,:)));
 title('单PN估计结果示意');
 subplot(1,2,2);
 plot(abs(channel_estimate_2(sim_num-5,:)),'r');
 title('双PN估计结果示意');
 figure;
 pn_freq = fft(channel_estimate_2(sim_num-5,:),Frame_len);
 plot(abs(pn_freq));
  
 SNR = estimate_SNR(recover_data(9:FFT_len+1:end-FFT_len),data_transfer(9:FFT_len+1:end-FFT_len))
 
 kk = 1;
 num = sim_num-14;
 for i = 9+(1:num)
     channel_off =  channel_estimate_1(i,:)- channel_estimate_2(i,:);
     channel_SNR(kk) = var(channel_estimate_2(i,:))/var(channel_off);
     kk = kk + 1;
 end
 chan_off_SNR = 10*log10(mean(channel_SNR))
 
