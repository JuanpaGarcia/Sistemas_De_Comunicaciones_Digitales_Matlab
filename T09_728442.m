clear all;
close all;
clc 

%Tarea 9

%%
%Efficient pulse design 


plots = 0;

Fs = 8000;
B = 1000;
Rs = 2000;
Rb = Rs;
Tp = 1/Rb; 
Ts = 1/Fs;
Entero_forzoso = Tp/Ts; %mp
energy = Tp;
%beta=(2*B/(Rs))-1; Beta=0
beta = 0;
D = 10;
type = 'rc';

pbase1 = rcpulse(beta,D,Tp,Ts,type,energy);

%Pulse Duration time D*Tp
%%
%Pulse 2
Fs = 8000;
Rs = 2000;
Entero_forzoso = Tp/Ts; %mp
beta = 0.2;
Rb = ((2*B)/Rs)-1;
Tp = 1/Rb; 
Ts = 1/Fs;
energy = Tp;
D = 10;
type = 'rc';

pbase2 = rcpulse(beta,D,Tp,Ts,type,energy);
%%
%Pulse 3
Fs = 4000;
Rs = 2000;
Entero_forzoso = Tp/Ts; %mp
beta = 0.8;
Rb = 2000;
Tp = 1/Rb; 
Ts = 1/Fs;
energy = Tp;
D = 6;
type = 'rc';
B=(Rb*(beta+1))/2;

pbase3 = rcpulse(beta,D,Tp,Ts,type,energy);
%%
%Exercise 2
if plots ==1
    wvtool(pbase1);
    wvtool(pbase2);
    wvtool(pbase3);
end

%%
%Exercise 3
beta = 0;
Fs = 1000;
Tp = 1/100;
D = 10;
Rs = 2000;
Ts = 1/Fs;
Entero_forzoso = Tp/Ts; %mp
mp = Tp/Ts; %mp
energy = Tp;
type = 'rc';

pbase = rcpulse(beta,D,Tp,Ts,type,energy);

V_bit = [1 0 1 1 0 0 1 1 1 1 1];

V_bit_polar = zeros(1,numel(V_bit)*mp);

counter = 0;
for i= 0 : numel(V_bit)-1
    if V_bit(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_bit_polar(counter*i+1) = value;
    counter = mp;
end
%%
% Polar NRZ LineCode halfsine
Polar_NRZ_sig = conv(pbase ,V_bit_polar);
plots = 1;
if plots == 1
    figure();
    plot(Polar_NRZ_sig);
    title('Transmitted signal');
end

%%
%Exercise 4
%Transmition 

 beta = 0;
 Fs = 96000;
 Ts = 1/Fs;
 Rs = 9600;
 D = 10;
 Rb = Rs;
 Tp = 1/Rb; 
 energy = sqrt(Tp);
 type = 'rc';
 
 pbase = rcpulse(beta,D,Tp,Ts,type,energy);





























