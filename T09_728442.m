clear all;
close all;
clc 

%Tarea 9

%%
%Efficient pulse design 

plots = 0;

Fs = 8000;
Rs = 2000;
Rb = Rs;
Tp = 1/Rb; 
Ts = 1/Fs;
Entero_forzoso = Tp/Ts; %mp
energy = Tp;
beta = 0;
D = 10;
type = 'rc';

pbase = rcpulse(beta,D,Tp,Ts,type,energy);

if plots == 1
    wvtool(pbase);
end

%Pulse Duration time D*Tp







