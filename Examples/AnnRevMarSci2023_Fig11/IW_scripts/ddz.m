function d=ddz(z)
% First derivative matrix for independent variable z. 
% 2nd order centered differences.
% Use one-sided derivatives at boundaries.
% z is assumed to be evenly spaced.
%
% Function provided here as an archive and downloaded on June 29 2022 from:
% https://blogs.oregonstate.edu/salty/matlab-tools-to-solve-the-viscous-taylor-goldstein-equation-for-both-instabilities-and-waves/
%
% Reference: 
% Smyth, W.D., J.N. Moum and J.D. Nash, 2011: 
%   ?Narrowband, high-frequency oscillations at the equator. Part II: 
%   Properties of shear instabilities?, J. Phys. Oceanogr. 41, 412-428.
%

% check for equal spacing
if abs(std(diff(z))/mean(diff(z)))>.000001
    disp(['ddz: values not evenly spaced!'])
    d=NaN;
    return
end

del=z(2)-z(1);N=length(z);

d=zeros(N,N);
for n=2:N-1
    d(n,n-1)=-1.;
    d(n,n+1)=1.;
end
d(1,1)=-3;d(1,2)=4;d(1,3)=-1.;
d(N,N)=3;d(N,N-1)=-4;d(N,N-2)=1;
d=d/(2*del);
return
end