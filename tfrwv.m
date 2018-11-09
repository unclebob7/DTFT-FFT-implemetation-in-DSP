function [tfr,t,f] = tfrwv(Sig,SampFreq,N)
%TFRWV	Wigner-Ville time-frequency distribution.
%	[TFR,T,F]=TFRWV(X,T,N,TRACE) computes the Wigner-Ville distribution
%	of a discrete-time signal X, 
%	or the cross Wigner-Ville representation between two signals. 
% 
%	X     : signal if auto-WV, or [X1,X2] if cross-WV.
%	T     : time instant(s)          (default : 1:length(X)).
%	N     : number of frequency bins (default : length(X)).
%	TRACE : if nonzero, the progression of the algorithm is shown
%	                                 (default : 0).
%	TFR   : time-frequency representation. When called without 
%	        output arguments, TFRWV runs TFRQVIEW.
%	F     : vector of normalized frequencies.
%
%	Example :
%	 sig=fmlin(128,0.1,0.4); tfrwv(sig);
% 
%	See also all the time-frequency representations listed in
%	the file CONTENTS (TFR*)

%	F. Auger, May-August 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

if (nargin == 0),
 error('At least one parameter required');
end;

Sig=hilbert(real(Sig));
SigLen=length(Sig);
if(nargin<3),
    N=512;
end
if(nargin<2),
    N=1;
end
if (N<0),
 error('N must be greater than zero');
end;
N=2^nextpow2(N);
N=min(N,SigLen);
tfr= zeros (N,SigLen);  
for icol=1:SigLen,
 taumax=min([icol-1,SigLen-icol,round(N/2)-1]);
 tau=-taumax:taumax; 
 indices= rem(N+tau,N)+1;
 tfr(indices,icol) =Sig(icol+tau) .* conj(Sig(icol-tau));
end; 
tfr= fft(tfr)/N; 
f=linspace(0,1,N)*SampFreq/2;
t=(0:SigLen-1)/SampFreq;
if (nargout==0)
set(gcf,'Position',[20 100 400 300]);
mesh(t,f,abs(tfr));
colorbar;
axis([min(t) max(t) min(f) max(f)]);
ylabel('f/Hz');
xlabel('t/s');
title('Wigner-Ville·Ö²¼')
end
end
