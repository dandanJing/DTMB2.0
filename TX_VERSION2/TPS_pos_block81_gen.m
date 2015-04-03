function [tps_pos] = TPS_pos_block81_gen(frame_num, tps_mode)
%tps_mode 0£º»¬¶¯40  1£º»¬¶¯27
block_len = 81;
if tps_mode == 0
    if mod(frame_num,2)==0
         start_pos = 41;
    elseif mod(frame_num,3)==1
         start_pos = 1;
    else
        start_pos = 81;
    end
elseif tps_mode == 1
    start_pos = 1+mod(frame_num-1,4)*27;
else
    start_pos = 1;
end
tps_pos = [start_pos:block_len:3888*8];