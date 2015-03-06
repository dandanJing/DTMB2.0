function h_out = channel_estimate_A(data_in, PN, FFT_length)

 %%op1 信道粗时域冲激响应
 fft_PN_R = fft(data_in, FFT_length);
 fft_PN = fft(PN, FFT_length);
 H_F =  fft_PN_R./fft_PN;
 
%平滑滤波
freq_thres = 3.6;
H_F(abs(fft_PN)<freq_thres)=0;
h_coarse = ifft(H_F);
mmse_weight = 0.99;
h_mmse_filter = channel_mmse_filter(h_coarse, mmse_weight);
h_out = h_mmse_filter.';
h_max = max(abs(h_out));
h_out(abs(h_out)<h_max*0.01)=0;
h_out(abs(h_out)~=0)=h_coarse(abs(h_out)~=0);

