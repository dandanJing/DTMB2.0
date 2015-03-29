%function: multipath channel model with ideal low pass filtering and M-times oversampling
%mp_mode: Medium-Echo Static (1--6); Violence-Echo Static (7--10); Brazil
%Static (11--15); dvb-t rayleigh/rician (16--17) 18£ºdvb-T2 0dB echo
%unit_delay: Ts=1/Fs
%M: oversampling factor
%updated by kewu peng on 2008-07-24
function h = multipath_new(mp_mode, unit_delay, M, dbg)

if (mp_mode>18 || mp_mode<1) h=1; return; end

path_scale(1:14,1:6)=0; %in db
path_delay(1:14,1:6)=0; %in us

path_scale(1,1:6)=[  -20    0  -20  -10  -14  -18];
path_delay(1,1:6)=[ -1.8    0 0.15  1.8  5.7   18];

path_scale(2,1:6)=[  -18    0  -20  -20  -10  -14];
path_delay(2,1:6)=[ -1.8    0 0.15  1.8  5.7   30];

path_scale(3,1:6)=[  -20    0  -14  -10  -20  -18];
path_delay(3,1:6)=[ -1.8    0 0.15  1.8  5.7   18];

% path_scale(4,1:6)=[    0  -10  -14  -18  -20  -20];
% path_delay(4,1:6)=[    0  0.2  1.9  3.9  8.2   15];
path_scale(4,1:6)=[    0  0  -80  -100  -150  -200];
path_delay(4,1:6)=[    0  50.5/7.56  1.9  3.9  8.2   15];
% path_scale(4,1:6)=[    0  -2  -3  -6  -8  -10];
% path_delay(4,1:6)=[    0  0.2  0.5  1.6  2.3   5];

path_scale(5,1:6)=[  -19    0  -22  -17  -22  -19];
path_delay(5,1:6)=[ -0.2    0 0.08 0.15  0.3  0.6];

path_scale(6,1:6)=[  -10  -20    0  -20  -10  -14];
path_delay(6,1:6)=[  -18 -1.8    0 0.15  1.8  5.7];

path_scale(7,1:6)=[  -20    0  -20  -10    0  -18];
path_delay(7,1:6)=[ -1.8    0 0.15  1.8  5.7   18];

path_scale(8,1:6)=[  -18    0  -20  -20  -10    0];
path_delay(8,1:6)=[ -1.8    0 0.15  1.8  5.7   30];
% path_scale(8,1:6)=[  -100    0  -100  -100  -100    0];
% path_delay(8,1:6)=[ -1.8    0 0.15  1.8  5.7   30];
% % path_delay(8,1:6)=[ -1.8    0 0.15  1.8  5.7   30.093];
% path_delay(8,1:6)=[ -1.8    0 0.15  1.8  5.7   30.16];

path_scale(9,1:6)=[ -5.1    0 -3.9 -3.8 -2.5 -1.3];
path_delay(9,1:6)=[ 0.07 0.52 0.60 0.85 2.75 3.23];

path_scale(10,1:6)=[    0 -0.5 -4.3 -4.4 -3.0 -1.8];
path_delay(10,1:6)=[ 0.43 0.52 0.85 1.37 2.75 3.23];

path_scale(11,1:6)=[    0 -13.8 -16.2 -14.9 -13.6 -16.4];
path_delay(11,1:6)=[    0  0.15  2.22  3.05  5.86  5.93];

path_scale(12,1:6)=[    0   -12    -4    -7   -15   -22];
path_delay(12,1:6)=[    0   0.3   3.5   4.4   9.5  12.7];

path_scale(13,1:6)=[ -2.8   0.0  -3.8  -0.1  -2.5  -1.3];
path_delay(13,1:6)=[  0.0 0.089 0.419 1.506 2.332 2.799];

path_scale(14,1:6)=[ -0.1  -3.8  -2.6  -1.3   0.0  -2.8];
path_delay(14,1:6)=[ 0.15  0.63  2.22  3.05  5.86  5.93];

path_scale(15,1:6)=[    0     0     0   -80   -80   -80];
path_delay(15,1:6)=[  0.0   1.0   2.0   3.0   4.0   5.0];

path_scale(18,1:2)=[    0     0  ];
path_delay(18,1:2)=[  0.0   51.5  ];

%tdb=[-7.8 -24.8 -15.0 -10.4 -11.7 -24.2 -16.5 -25.8 -14.7 -7.9 -10.6  -9.1 -11.6 -12.9 -15.3 -16.5 -12.4 -18.7 -13.2 -11.7]; t=10.^(tdb/20); norm(t)
%attenuation, delay (in us) and phase (in rad) values
x=[...
 1 0.057662 1.003019 4.855121;
 2 0.176809 5.422091 3.419109;
 3 0.407163 0.518650 5.864470;
 4 0.303585 2.751772 2.215894;
 5 0.258782 0.602895 3.758058;
 6 0.061831 1.016585 5.430202;
 7 0.150340 0.143556 3.952093;
 8 0.051534 0.153832 1.093586;
 9 0.185074 3.324866 5.775198;
10 0.400967 1.935570 0.154459;
11 0.295723 0.429948 5.928383;
12 0.350825 3.228872 3.053023;
13 0.262909 0.848831 0.628578;
14 0.225894 0.073883 2.128544;
15 0.170996 0.203952 1.099463;
16 0.149723 0.194207 3.462951;
17 0.240140 0.924450 3.664773;
18 0.116587 1.381320 2.833799;
19 0.221155 0.640512 3.334290;
20 0.259730 1.368671 0.393889];

%--------%get the transfer function of multipath channel
%system impulse response for ideal low pass filter is sin(pi*t/Ts)/(pi*t/Ts)=sinc(t/Ts)
Ts=unit_delay; T=Ts/M;
if (mp_mode==16 || mp_mode==17)
  scale=x(:,2)'.*exp(-j*x(:,4)'); %rou=rou/norm(rou);
  delay=x(:,3)';
else
  delay=path_delay(mp_mode,:);
  scale=path_scale(mp_mode,:); scale=10.^(scale/20);
end

if (mp_mode==17) %dvb-t rician channel
  %rou(0)^2=sum(rou(i)^2)*10 (K=10 for rician factor)
  scale=[scale sqrt(10)*norm(scale)];
  delay=[delay 0];
end

m1=round(min(delay)/T); m2=round(max(delay)/T);
% N1=m1-4*M; N2=m2+4*M;
N1=m1- 21*M; N2=m2+21*M;

h=zeros(1,N2-N1+1);
for (i1=1:length(scale))
  t=(N1:N2)*T-delay(i1);
  h=h+scale(i1)*sinc(t/Ts);
end
h=h/norm(h);

if (dbg) figure; subplot(2,1,1); plot(real(h),'b.-'); subplot(2,1,2); plot(abs(fft(h,4096)),'b.-'); axis tight; end
return