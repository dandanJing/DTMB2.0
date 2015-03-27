function tps_pos = scatter_tps_position_gen(frame_num, tps_mode)
block_len = 96;
start_pos = 1+mod(frame_num-1,4)*24;
%start_pos = block_len - mod(frame_num-1,4)*24;
tps_pos = start_pos:block_len:3888*8;