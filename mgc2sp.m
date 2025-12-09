function lnH=mgc2sp(mgcep, alfa, gamma, fftsize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mgc2sp.m
% Created on: Jan 05, 2010
% Author: Fabio Tesser
% Institution: CNR-ISTC, Padova - Italy
% Email: fabio.tesser@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MGC2SP calculates the log (ln) magnitude spectrum from mel-generalized cepstral coefficients mgcep
% mgc2sp(mgcep, alfa, gamma, fftsize)
%
% fftsize = N
%   N even:	(N+2)/2 points are returned with
% 			the first and last being real
%   N odd:	(N+1)/2 points are returned with the
% 			first being real

if gamma ~= 0
    error('gamma ~= 0')
    exit
end

num_samples=size(mgcep,1);
mgspec_order=size(mgcep,2);

for sample=1:num_samples
    % not necessary to compute all (real simmetric)
    % 1+fix(n/2)
    K= 1+fix(fftsize/2);
    for k=0:K-1 
        tmp=0;
        % not necessary to compute all %n=0:fftsize-1 if fftsize > mgspec_order
        % but only mgspec_order
        for n=0:mgspec_order-1 
            z_inv=exp(((-2*pi*i)/fftsize)*k);
            warped_z_inv = freq_trans(z_inv, alfa);
            tmp=tmp+(mgcep(sample,n+1)*(warped_z_inv^n));
        end
        lnH(sample,k+1)=tmp;
    end
end

    

