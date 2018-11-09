%%读取数据
clear;
clc;
fid=fopen('Force.Z.SAC','r','ieee-le');
A=fread(fid,[70,1],'float32');
B=fread(fid,[40,1],'int32');
C=char(fread(fid,[1,192],'char'));
HR=fread(fid,'float32');
A(A==-12345.0)=NaN;
B(B==-12345)=NaN;
fclose(fid);

xn=HR';
xn=[xn zeros(1,27)];
xn=xn';
N=length(xn);
fs = 2;   
t = 0:N-1/fs;
fre_bins=1024;
threshold=10;
threshold_1=3;
tic;  %计算时间
[TF_st,TF_wvd,tfr,f,t]=gyh_stft_wvd(xn,t,fs,fre_bins,threshold,threshold_1);
toc;
figure(1)
imagesc(t,f,abs(tfr));
set(gca,'YDir','normal');
axis([0 3000 0 0.4]);
colorbar
colormap ('Jet');
ylabel('Frequency/Hz');
xlabel('Time/Sec');
title('原始信号的 stft-wvd 分布');

