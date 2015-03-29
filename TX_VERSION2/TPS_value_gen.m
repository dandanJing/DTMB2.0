m_init = [ 1 1 1 1 1 1 1 1 1 1 1];
m_len = 2^11 - 1;
m_sequence = zeros(1,m_len);
for kk = 1: m_len
    temp = mod(m_init(9)+m_init(11),2);
    m_sequence(kk) = m_init(end);
    m_init = [temp m_init(1:end-1)];
end
tps_value_temp = 1- 2*m_sequence(1:48*8);
tps_value = (tps_value_temp+j*tps_value_temp)*1024*3.5*10^2;
save('TPS_symbol','tps_value');