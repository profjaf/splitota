function [x,bw] = vkSplitOneOrd(y,theta,fs,r,filtord)
%
% Type: [x,bw] = vkSplitOneOrd(y,f,fs,r,filtord)
%
% Inputs:
%
% y    := N x 1, data vector
% theta    := N x 1, angular position vector [rad] 
% fs   := 1 x 1, sampling frequency [Hz]
% r    := 1 x 1, weighting factor for the structural equation (typical value in 100's or 1000's)
% filtord := 1 x 1, filter order, one less than the number of poles.  Can take values 1 or 2.
%
% Outputs:
%
% x    := 2N x 1, output order vector minimizing the objective function
%         
% bw   := Nord x 1, -3 dB filter bandwidth for each filter [Hz]
%
% This function is for extracting a single order at a time (not simultaneous orders).
%Cite as:
% Improved Multi-order Vold-Kalman Filter for Order Tracking Analysis Using Split Cosine and Sine Terms
% Written by jaf, July 2023.
%
if ~isvector(y)
 error('y must be a vector.');
end
y=y(:); N = length(y);
if ~isvector(theta)
 error('theta must be a vector.');
end
N2 = length(theta);
if N2 ~= N
 error('length(f) must = length(y).');
end
if ~ismember(filtord,[0,1,2])
 error('filtord must be element of [0,1,2].');
end
%
dt = 1/fs;
% this is part of the code only
% full code available upon request

% [R,p] = chol(GA);%in fact GA is symmetric positive definite, so it will
% be eventually solved by Cholesky factorization
% p is zero when GA is symmetric positive definite
x = (GA\yy); %this gives the required solution
%x = lsqminnorm(GA,yy); %failed to give good results
%x = 1*lsqr(GA,yy); %overfit data equations but obtain sinusoidal envelope



%
%%