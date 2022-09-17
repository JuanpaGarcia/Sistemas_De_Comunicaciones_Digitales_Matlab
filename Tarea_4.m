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
x = conv(s,pbase);  % Verifica las dimensiones de s,p y x. ¿Qué relación hay entre éstas? 
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
bits=b(:);   % Vector de bits concatenado 
%%
%Recuperar la imagen

recuperado = zeros(33,33,'uint32'); %crear espacio apra imagen
%ciclo para cargar y convertir valores en la matriz nuewva
counter = 8; %variable para avanzar en el vector de bits
for i = 1 : 33
    for j = 1: 33
        recuperado(j,i) = bi2de(b(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end
imshow(uint8(recuperado)) %imagen recuperada

%%
%Generar los pulsos 

mp = 20;

V_16bit = b(1:16);

%definir pulsos
%unipolar NRZ
UPNRZ = ones(1,mp);
stem(UPNRZ)
wvtool(UPNRZ) 

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
V_16bit_polar(1:mp:end) = V_16bit*2-1; 


V_16bit_Unipolar = zeros(1,numel(V_16bit)*mp);
V_16bit_Unipolar(1:mp:end) = V_16bit;
%%
%Señal con Unipolar NRZ
Unipolar_NRZ_Sig = conv(UPNRZ ,V_16bit_Unipolar);
plot(Unipolar_NRZ_Sig)

%%
%Selak Polar NRZ
Polar_NRZ_Sig = conv(UPNRZ ,V_16bit_polar);
plot(Polar_NRZ_Sig)

%%
%Señal Polar RZ
Polar_RZ_sig = conv(PRZ, V_16bit_polar);
plot(Polar_RZ_sig)

%%
%Señal Bipolar NRZ
Bipolar_NRZ = conv(UPNRZ, V_16bit_polar);
plot(Bipolar_NRZ);

%%
%Señal Manchester
Manchester = conv(Manchester, V_16bit_polar);
plot(Manchester);







