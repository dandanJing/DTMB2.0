function chlen =chan_len_estimate_forpower(h_in)
h_power = norm(h_in);
thresh = sqrt(0.95);
h_max = max(abs(h_in));
h_pos = find(abs(h_in)>0.05*h_max);
for k = h_pos(end):length(h_in)
    h_power_temp = norm(h_in(1:k))/h_power;
    if h_power_temp > thresh
        break;
    end
end
chlen = k;