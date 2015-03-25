function [tps_pos tps_value] = TPS_gen(frame_num, tps_mode)
m_init = [ 1 1 1 1 1 1 1 1 1 1 1];
m_len = 2^11 - 1;
m_sequence = zeros(1,m_len);
for kk = 1: m_len
    temp = mod(m_init(9)+m_init(11),2);
    m_sequence(kk) = m_init(end);
    m_init = [temp m_init(1:end-1)];
end
tps_value_temp = 1- 2*m_sequence(1:48*8);

block_len = 96;
start_pos = block_len - mod(frame_num-1,4)*24;
scater_pos = start_pos:block_len:3888*8;
conti_len = 30;
scater_min = min(scater_pos);
scater_max = max(scater_pos);
conti_pos = [1:conti_len 3888*8-conti_len+1:3888*8];
if(find(conti_pos == scater_min))
    conti_pos(conti_pos == scater_min) = conti_len+1;
elseif(find(conti_pos == scater_max))
    conti_pos(conti_pos == scater_max) = 3888*8-conti_len;
end
tps_pos = [scater_pos conti_pos];
tps_pos = sort(tps_pos);
tps_value = (tps_value_temp+j*tps_value_temp)*1024;