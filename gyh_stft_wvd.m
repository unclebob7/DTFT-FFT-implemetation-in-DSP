function [TF_st,TF_wvd,tfr,f,T]=gyh_stft_wvd(x,t,fs,fre_bins,threshold,threshold_1)

%   x     : signal. ��������ʵֵ
% 	t     : time instant(s) (1:length(x))         .
% 	fre_bins     : number of frequency bins .(2^n)
%   fs��   �����ʣ�
%   threshold    ��ֵ���������������
%   threshold_1  ��ֵ1:���������ź����Ͻ���� ��1 ~ 2��
%   F 0~0.5
Nw=2*fre_bins;
h=window(@gausswin,Nw-1);

[tfr_st,T,f] = tfrstft(x,t,Nw,h);

TF_st=abs(tfr_st);
TF_st=TF_st(1:fre_bins,:);


[ tfr_wvd,T,f] = tfrwv( x,fs,fre_bins);
TF_wvd=abs(tfr_wvd);

%F=0:0.5/(fre_bins-1):0.5;
%% STFT�����һ������
ST=TF_st;
max_stft=max(max(ST));
ST=ST/max_stft;%��һ��
[i,j,v]=find(ST>=1);
[a,b,c]=find(ST);%�ҳ�ST������������C
ST(ST==0)=min(c);%��ST������Ϊ�����ȫ���滻�ɷ������Сֵ
%% WVD�����һ������
tfr=TF_wvd;
max_wvd=tfr(i,j);
A=tfr/max_wvd;
%-----------------------------------------------%
B=A./ST;
tfr(B>threshold)=0; %����������
B(B>threshold)=1;
B(B<threshold_1)=1;
tfr=tfr./B;
%% ��Сֵ��
% TF_wv=TF_wvd;
% B=TF_wv./TF_st;
% TF_wv(B>1)=TF_st(B>1);
% tfr=TF_wv;
%% �ݵ���ϵ����
%tfr=(TF_wvd.^0.5).*(TF_st.^0.5);
%% ��ֵ����
 TF_st(TF_st<2)=0;
 TF_st(TF_st>1)=1;
 tfr=TF_wvd.*TF_st;
