clear all;
close all;
clc 
%%
mp = 11; % muestras por pulso 
pbase = triang(mp); % ejemplo con pulso base triangular de mp muestras
b = [1 0 1 0 1 1 0 0 1 0 1 0 1 1 0 0]; 
s = zeros(1,numel(b)*mp);
s(1:mp:end) = b *2-1;
stem(s)
%%
x = conv(s,pbase);  % Verify dimentions
plot(x) 
hold on 
stem(x,'LineStyle','none') 

%%
%Ejercicio1
%load lena512.mat
%tuve que cargar la imagen manual porque no es soportada en mi version de
%matlab
lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(252:284,318:350); 
imshow(uint8(lenarec)) 

%%
b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector 
%%
%Image recovery

recuperado = zeros(33,33,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : 33
    for j = 1: 33
        recuperado(j,i) = bi2de(b(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %recovered image

%%
%Pulse generation

mp = 20;

V_16bit = b(1:16);

%Pulses
%unipolar NRZ
UPNRZ = ones(1,mp);
stem(UPNRZ)
wvtool(UPNRZ) 
%%
%Polar RZ
PRZ= zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        PRZ(i) = 1;  
    else
        PRZ(i) = 0; 
    end
end
stem(PRZ)
wvtool(PRZ) 

%%
%Manchester
Manchester = zeros(1,mp); 
for i = 1:mp
    if i <= mp/2
        Manchester(i) = 1;  
    else
        Manchester(i) = -1; 
    end
end
stem(Manchester)
wvtool(Manchester) 




%%
%Generación de las señales

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

V_16bit_Unipolar = zeros(1,numel(V_16bit)*mp);
V_16bit_Unipolar(1:mp:end) = V_16bit;
%%
Fs = 9600;
%Señal con Unipolar NRZ
Unipolar_NRZ_Sig = conv(UPNRZ ,V_16bit_Unipolar);
plot(Unipolar_NRZ_Sig)
title('Signal on Unipolar NRZ line code');
figure;
pwelch(Unipolar_NRZ_Sig,Fs,'power');

%%
%Selak Polar NRZ
Polar_NRZ_Sig = conv(UPNRZ ,V_16bit_polar);
plot(Polar_NRZ_Sig)
title('Signal on Polar NRZ line code');

%%
%Señal Polar RZ
Polar_RZ_sig = conv(PRZ, V_16bit_polar);
plot(Polar_RZ_sig)
title('Signal on Polar RZ line code');

%%
%Señal Bipolar NRZ
Bipolar_NRZ = conv(UPNRZ, V_16bit_polar);
plot(Bipolar_NRZ);
title('Signal on Bipolar NRZ line code');

%%
%Señal Manchester
Manchester = conv(Manchester, V_16bit_polar);
plot(Manchester);
title('Signal on manchester line code');






