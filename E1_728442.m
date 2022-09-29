clear all;
close all;
clc 

%Examen 1

Tp = 1/(9420);  %Periodo pulso
Fs = 2*(1/Tp);  %Nyquist teorema 
Ts = 1/Fs;
P = [0.0 0.242 0.442 0.542 0.742 0.842 0.942 0.99 0.99 0.942 0.842 0.742 0.542 0.442 0.242];

%%
wvtool(P);
pwelch(P,[],[],[],Fs,'power');


%%
%Parte 6
lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(221:284, 287:350); 
imshow(uint8(lenarec)) 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector 

%Capture first 16 bits

V_16bit = b(1:16);

mp = 20;

Polar = ones(1,mp);

V_16bit_polar = zeros(1,numel(V_16bit)*mp);
%V_16bit_polar(1:mp:end) = V_16bit*2-1; 
counter = 0;
for i= 0 : numel(V_16bit)-1
    if V_16bit(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_16bit_polar(counter*i+1) = value;
    counter =20;
end

Polar_NRZ_Sig = conv(P ,V_16bit_polar);
plot(Polar_NRZ_Sig)
title('Codigo de linea Polar NRZ');
%%
pwelch(Polar_NRZ_Sig,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Spectral power density');


%%
%filtro 
F_filtro_max = Fs/2;
Fc_ejer = 1987.03/Fs *pi; %Fc = 1987 hz
Fcorte = Fc_ejer + 0.021; %edad 21
Fcorte = Fcorte /pi;

f = [0 Fcorte Fcorte 1];
m = [1 1 0 0]; o = 60;
LPF = fir2(o,f,m);
fvtool(LPF);  % Espectro en frecuencia del filtro.

%%
%Parte 9
%Capture first 16 bits

V_200bit = b(1:200);

mp = length(P);

V_200bit_polar = zeros(1,numel(V_200bit)*mp);

counter = 0;
for i= 0 : numel(V_200bit_polar)-1
    if V_200bit_polar(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_200bit_polar(counter*i+1) = value;
    counter =mp;
end

Polar_NRZ_200Sig = conv(P ,V_200bit_polar);
plot(Polar_NRZ_200Sig)
title('Codigo de linea Polar NRZ');
%%

Salida_filtro = conv(Polar_NRZ_200Sig ,LPF);
%%
plot(Salida_filtro);
%%
pwelch(Salida_filtro,[],[],[],Fs,'power');  % PSD of Unipolar NRZ
title('Spectral power density');
