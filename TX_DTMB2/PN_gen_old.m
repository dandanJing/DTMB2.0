function PN=PN_gen_old;
PN_Rom_Create_420;
PN_rom = 1 - 2*PN_rom;
PN_total=432;
PN_len=255;
PN_cyclic_Len=PN_total-PN_len;
frame_num_md =  1;
PN_Selected = PN_rom(frame_num_md, 2:PN_len+1);
PN = [PN_Selected(PN_len - PN_cyclic_Len + 1:end),PN_Selected];
PN=PN*1024;
PN=PN+j*PN;