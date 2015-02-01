function [sc_h1, sc_h2, sc_ha] = h_estimate_A(h_pre2, h_pre1, pos, thresh)

h_pre_max = max(abs(h_pre1));
abs_h_pre = abs(h_pre1);
alpha = 0.97;

if(pos>=4)
	fc_dif = h_pre1 - h_pre2;
	fc_dif( abs_h_pre<h_pre_max*thresh) = 0;
else
	fc_dif = zeros(1, 1);
end

sc_h1 = h_pre1 + alpha*fc_dif;
sc_h2 = h_pre1 + 2*alpha*fc_dif;
sc_ha = h_pre1 + 1.5*alpha*fc_dif;