%%DTMB2.0数据接收 帧头432，帧体3888*8，TPS 48*8,64QAM
clear all,close all,clc

fid = fopen('dtmb2_data.txt', 'r');
Send_data_temp = fscanf(fid, '%03x');
fclose(fid);


