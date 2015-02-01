function h_out = channel_estimate(data_in, PN)

 %%op1 信道粗时域冲激响应
 fft_PN_R = fft(data_in,2048);
 fft_PN = fft(PN,2048);
 H_F =  fft_PN_R./fft_PN;
      
%平滑滤波
freq_thres = 3.6;
mmse_weight = 0.99;
H_F(abs(fft_PN)<freq_thres)=0;
h_coarse = ifft(H_F);
h_mmse_filter = channel_mmse_filter(h_coarse, mmse_weight);

h_out = h_mmse_filter.';


