function PN=PN_gen(frame_num, pn_mode);
if pn_mode == 0 %%PN432
    PN_Rom_Create_420;
    PN_rom = 1 - 2*PN_rom;
    PN_total=432;
    PN_len=255;
    PN_cyclic_Len=PN_total-PN_len;
    frame_num_md =  mod(frame_num-1,225)+1;
    PN_Selected = PN_rom(frame_num_md, 2:PN_len+1);
    PN = [PN_Selected(PN_len - PN_cyclic_Len + 1:end),PN_Selected];
else
    PN_Rom_Create_945;
	PN_rom = 1 - 2*PN_rom;
    frame_num_md =  mod(frame_num-1,200)+1;
	PN = PN_rom(frame_num_md, :);
end
PN=PN*1024*2;
PN=PN+j*PN;