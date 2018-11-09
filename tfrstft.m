function [tfr,t,f] = tfrstft(x,t,N,h,trace);

% TFRSTFT Short time Fourier transform.

% 	[TFR,T,F]=TFRSTFT(X,T,N,H,TRACE) computes the short-time Fourier

% 	transform of a discrete-time signal X. 

% 
% 	X     : signal.

% 	T     : time instant(s)          (default : 1:length(X)).

% 	N     : number of frequency bins (default : length(X)).

% 	H     : frequency smoothing window, H being normalized so as to

% 	        be  of unit energy.      (default : Hamming(N/4)). 

% 	TRACE : if nonzero, the progression of the algorithm is shown

% 	                                 (default : 0).

% 	TFR   : time-frequency decomposition (complex values). 

%The frequency axis is graduated from -0.5 to 0.5.

% 	F     : vector of normalized frequencies.

% 
% 	Example :

% 	 sig=[fmconst(128,0.2);fmconst(128,0.4)]; tfr=tfrstft(sig);

% 	 subplot(211); imagesc(abs(tfr));

% 	 subplot(212); imagesc(angle(tfr));

% 
% 	See also all the time-frequency representations listed in

% 	 the file CONTENTS (TFR*)

% 
% 	F. Auger, May-August 1994, July 1995.

% 	Copyright (c) 1996 by CNRS (France).

% 
% 	------------------- CONFIDENTIAL PROGRAM -------------------- 

% 	This program can not be used without the authorization of its

% 	author(s). For any comment or bug report, please send e-mail to 

% 	f.auger@ieee.org 

%% TFRSTFT��ʱ����Ҷ�任��

% [TFR��T��F]= TFRSTFT��X��T��N��H��TRACE�������ʱ����Ҷ

% ��ɢʱ���ź�X�ı任

% 
% X���źš�

% T��ʱ�̣�S����Ĭ��ֵ��1�����ȣ�X������

% N��Ƶ�ʶ�����Ĭ��ֵ�����ȣ�X������

% H��Ƶ��ƽ�����ڣ�H���飬��
% ��������λ�� ��Ĭ��ֵ��������N/4������

% TRACE��������㣬���㷨�Ľ�չ��ʾ
% ��Ĭ��ֵ��0����

% TFR��ʱƵ�ֽ⣨����ֵ����
%��Ƶ����-0.5��ҵ��0.5��

% F����׼��Ƶ�ʵ�������

% 
% ���磺

% SIG=[fmconst��128,0.2��; fmconst��128,0.4��]; TFR= tfrstft��SIG��;

% ������211��;��imagesc��ABS��TFR����;

% ������212��;��imagesc���ǣ�TFR����;

% 
% �����������г���ʱ�� - Ƶ�ʱ�ʾ

% ���ļ������ݣ�TFR*��

% 
% F.��Ъ��5��1994��8�£�1995��7�¡�

% ��Ȩ���У�c��1996���ɷ������ҿ�ѧ�о����ģ���������

% 
% -------------------����PROGRAM--------------------
% ���������δ����Ȩʹ����
% ���ߣ�S���������κ��������󱨸棬�뷢�͵����ʼ���
% f.auger@ieee.org
%%

%wait=waitbar(0,'STFT��ֵ������! ���Ժ�...');


[xrow,xcol] = size(x);
if (nargin < 1),
 error('At least 1 parameter is required');

elseif (nargin <= 2),
 N=xrow;
end;

hlength=floor(N/4);

hlength=hlength+1-rem(hlength,2);


if (nargin == 1),
 t=1:xrow; 
h = window(hlength); 
trace=0;

elseif (nargin == 2) || (nargin == 3),
 h = window(hlength); 
trace = 0;
elseif (nargin == 4),
 trace = 0;

end;


if (N<0),
 error('N must be greater than zero');

end;

[trow,tcol] = size(t);

if (xcol~=1),
 error('X must have one column');

elseif (trow~=1),
 error('T must only have one row'); 

elseif (2^nextpow2(N)~=N),
 fprintf('For a faster computation, N should be a power of two\n');

end; 


[hrow,hcol]=size(h); 
Lh=(hrow-1)/2; 

if (hcol~=1)||(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');

end;


h=h/norm(h);


tfr=zeros(round(N/2),tcol) ; 
tfr1=zeros(N,tcol);
if trace, disp('Short-time Fourier transform');
 end;

for icol=1:tcol,
    ti=t(icol);
tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);

 indices= rem(N+tau,N)+1; 

 if trace, disprog(icol,tcol,10); 
end;
 
tfr1(indices,icol)=x(ti+tau,1).*conj(h(Lh+1+tau));
 
%waitbar(icol/tcol,wait);

end;

tfr1=fft(tfr1); 
tfr=tfr1(1:round(N/2),:);
if trace, fprintf('\n'); 
end;


if (nargout==0),
 tfrqview(abs(tfr).^2,x,t,'tfrstft',h);

elseif (nargout==3),
 f=fftshift(0.5*(-N:N-1)/N)';

% close(wait);


end;
