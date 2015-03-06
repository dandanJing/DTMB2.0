function h_out = channel_estimate(data_in, PN, FFT_length,min_zeors_thresh)
%min_zeors_thresh 幅度小于此百分比阈值的径置零

 %%op1 信道粗时域冲激响应
 fft_PN_R = fft(data_in, FFT_length);
 fft_PN = fft(PN, FFT_length);
 H_F =  fft_PN_R./fft_PN;
      
%平滑滤波
freq_thres = 3.6;
mmse_weight = 0.99;

H_F(abs(fft_PN)<freq_thres)=0;
h_coarse = ifft(H_F);
h_mmse_filter = channel_mmse_filter(h_coarse, mmse_weight);
h_max = max(abs(h_mmse_filter));
h_out = h_mmse_filter.';
h_out(abs(h_out)<h_max*min_zeors_thresh )=0;
h_out(length(PN)+1:end)=0;


