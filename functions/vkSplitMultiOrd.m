function [x,bw] = vkSplitMultiOrd(y,theta,fs,r,filtord)
% 
% 
% Inputs:
%
% y       := N x 1, data vector
% theta   := N x Nord, matrix of radian position vectors for each order in columns [rad]
%            (y and f have the same number of rows)
% fs      := 1 x 1, sampling frequency [Hz]
% r       := Nord x 1, vector of weighting factors for the structural equations
%            (typical value in 100's or 1000's)
% filtord := 1 x 1, filter order, one less than the number of poles.  Can take values 1 or 2.
%
% Outputs:
%
% x    := N x Nord, output order vectors (in columns)
%         minimizing the multi-order objective function
%         x is the zero-peak complex amplitude envelope
% bw   := Nord x 1, -3 dB filter bandwidth for each filter [Hz]
%Cite as:
% Improved Multi-order Vold-Kalman Filter for Order Tracking Analysis Using Split Cosine and Sine Terms
% Written by jaf, July 2023.
%
%
if ~isvector(y)
 error('y must be a vector.');
end
y=y(:); N = length(y);
[N2,Nord] = size(theta);
if N2 ~= N
 error('size(f,1) must = length(y).');
end
if ~ismember(filtord,[0,1,2])
 error('filtord must be element of [1,2].');
end
if ~isvector(r)
 error('r must be a vector.');
end
r=r(:);
Nord2=length(r);
if Nord2 ~= Nord
 error('length(r) must = size(f,2).');
end
%
dt = 1/fs;
% this is part of the code only
% full code available upon request
%
% theta=2*pi*cumsum(f,1)*dt;%multi column angular position
% theta=theta-theta(1,:);%subtract the first angle since the angle must be zero at the begining

tmp = H\cy;%using direct solution, it is more accurate but not always guranteed
           %due to ill-conditioned system especially for large number of
           %points, where computer accuracy may be insufficient to provide
           %the required solution
%[tmp,fl,rr,it,rv] = pcg1(H,cy,tol,maxit,M1,x0); % external m-function
x = reshape(tmp,2*N,Nord);
