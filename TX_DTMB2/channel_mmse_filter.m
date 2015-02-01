function sc_h_out = channel_mmse_filter(sc_h_in, alpha)

fft_size = length(sc_h_in);
t_max = 0;
sc_h_out = complex(zeros(fft_size,1), zeros(fft_size,1));

t_max1 = max(abs(real(sc_h_in)));
t_max2 = max(abs(imag(sc_h_in)));
t_max = max(t_max1, t_max2);

y_c = (1-alpha)/alpha*t_max^2;

for i=1:fft_size
	t_I = real(sc_h_in(i))^2;
	t_Q = imag(sc_h_in(i))^2;
	
	yi = t_I/(t_I+y_c)*real(sc_h_in(i))/alpha;
	yq = t_Q/(t_Q+y_c)*imag(sc_h_in(i))/alpha;
    sc_h_out(i) = complex(yi, yq);
end
