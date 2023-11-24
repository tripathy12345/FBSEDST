%%%% Code for Fourier Bessel domain Discrete Stockwell transform%%%%
%if anyone is interested to use FBSE-DST, he/she can cite the following%%
%%%paper
%%Dash, Shaswati, Samit Kumar Ghosh, Rajesh Kumar Tripathy, Ganapati Panda, and Ram Bilas Pachori.
%"Fourier-Bessel domain based discrete Stockwell transform for the analysis of non-stationary signals."
%In 2022 IEEE India Council International Subsections Conference (INDISCON), pp. 1-6. IEEE, 2022.

clc;
clear all;
close all;
load sig.mat;
fs=256;
x=sig';
[a3, alfa, fre, y]=FBSE_S_transform(x,fs);
t=1:length(x);
contour(t,fre,y,10)
