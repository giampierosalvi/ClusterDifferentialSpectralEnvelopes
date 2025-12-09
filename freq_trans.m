%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% freq_trans.m
% Created on: Jan 05, 2010
% Author: Fabio Tesser
% Institution: CNR-ISTC, Padova - Italy
% Email: fabio.tesser@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FREQ_TRANS calculates the warped frequency Z using an all pass biliniar transformation
% freq_trans(zeta_inv,alfa)
% warp frequency inv zeta 

function warped_zeta_inv=freq_trans(zeta_inv,alfa)
warped_zeta_inv=(zeta_inv-alfa)/(1-alfa*zeta_inv);
