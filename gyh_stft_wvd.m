function [TF_st,TF_wvd,tfr,f,T]=gyh_stft_wvd(x,t,fs,fre_bins,threshold,threshold_1)

%   x     : signal. 列向量、实值
% 	t     : time instant(s) (1:length(x))         .
% 	fre_bins     : number of frequency bins .(2^n)
%   fs：   采样率；
%   threshold    阀值：消大多数交叉项
%   threshold_1  阈值1:用来消除信号项上交叉项； （1 ~ 2）
%   F 0~0.5
Nw=2*fre_bins;
h=window(@gausswin,Nw-1);

[tfr_st,T,f] = tfrstft(x,t,Nw,h);

TF_st=abs(tfr_st);
TF_st=TF_st(1:fre_bins,:);


[ tfr_wvd,T,f] = tfrwv( x,fs,fre_bins);
TF_wvd=abs(tfr_wvd);

%F=0:0.5/(fre_bins-1):0.5;
%% STFT数组归一化过程
ST=TF_st;
max_stft=max(max(ST));
ST=ST/max_stft;%归一化
[i,j,v]=find(ST>=1);
[a,b,c]=find(ST);%找出ST数组里非零的数C
ST(ST==0)=min(c);%将ST数组里为零的数全部替换成非零的最小值
%% WVD数组归一化过程
tfr=TF_wvd;
max_wvd=tfr(i,j);
A=tfr/max_wvd;
%-----------------------------------------------%
B=A./ST;
tfr(B>threshold)=0; %消除交叉项
B(B>threshold)=1;
B(B<threshold_1)=1;
tfr=tfr./B;
%% 最小值法
% TF_wv=TF_wvd;
% B=TF_wv./TF_st;
% TF_wv(B>1)=TF_st(B>1);
% tfr=TF_wv;
%% 幂调节系数法
%tfr=(TF_wvd.^0.5).*(TF_st.^0.5);
%% 二值化法
 TF_st(TF_st<2)=0;
 TF_st(TF_st>1)=1;
 tfr=TF_wvd.*TF_st;
