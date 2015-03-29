function [tps_pos] = TPS_pos_gen(frame_num, tps_mode)
%tps_mode 0:total_pos 1:scatter_pos 2:continuel_pos
block_len = 96;
start_pos = 1+mod(frame_num-1,4)*24;
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

if tps_mode == 1
    tps_pos = scater_pos;
elseif tps_mode == 2
    tps_pos = conti_pos;
else
    tps_pos = [scater_pos conti_pos];
    tps_pos = sort(tps_pos);
end