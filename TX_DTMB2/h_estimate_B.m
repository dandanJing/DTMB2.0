function [sc_h2, sc_ha] =  h_estimate_B(sc_h0, sc_h1, pos, thr_h_pred)

t_max = max(abs(sc_h1));
abs_h1 = abs(sc_h1);
alpha = 0.97;

if(pos>=4)
	fc_dif = sc_h1 - sc_h0;
	fc_dif(abs_h1<t_max*thr_h_pred) = 0;
else
	fc_dif = zeros(1, 1);
end

sc_h2 = sc_h1 + alpha*fc_dif;
sc_ha = sc_h1 + 0.5*alpha*fc_dif;
