clear all;
close all;
clc 

%Tarea 2    

[S,Fs] = audioread('spring_HIFI.wav');
seconds = 10;

for i = 1: (Fs*seconds)
    y(i) = S(i);
end
figure();
spectrogram(y);
title('Espectrograma de canción original');

%%
k=4;                                 % Cantidad de bits para Cuantizar
% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2
xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int
decVal = binaryVectorToDecimal(xq_bin,'MSBFirst');
decVal_song = int16(decVal);

filename = 'Audio_4b.wav';
audiowrite(filename,decVal_song,Fs,'BitsPerSample',8,...
'Comment','This is my new audio file.');
figure();
spectrogram(decVal);
title('Espectrograma de 4Bits');

%%
k=6;                                 % Cantidad de bits para Cuantizar
% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2
xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int
decVal = binaryVectorToDecimal(xq_bin,'MSBFirst');
decVal_song = uint8(decVal);

filename = 'Audio_6b.wav';
audiowrite(filename,decVal_song,Fs,'BitsPerSample',8,...
'Comment','This is my new audio file.');
figure();
spectrogram(decVal);
title('Espectrograma de 6Bits');

%%
k=8;                                 % Cantidad de bits para Cuantizar
% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2
xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int
decVal = binaryVectorToDecimal(xq_bin,'MSBFirst');
decVal_song = uint8(decVal);

filename = 'Audio_8b.wav';
audiowrite(filename,decVal_song,Fs,'BitsPerSample',8,...
'Comment','This is my new audio file.');
figure();
spectrogram(decVal);
title('Espectrograma de 8Bits');


%%
k=10;                                 % Cantidad de bits para Cuantizar
% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2
xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int
decVal = binaryVectorToDecimal(xq_bin,'MSBFirst');
decVal_song = int16(decVal);

filename = 'Audio_10b.wav';
audiowrite(filename,decVal_song,Fs,'BitsPerSample',16,...
'Comment','This is my new audio file.');
figure();
spectrogram(decVal);
title('Espectrograma de 10Bits');


%%
k=12;                                 % Cantidad de bits para Cuantizar
% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2
xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int
decVal = binaryVectorToDecimal(xq_bin,'MSBFirst');
decVal_song = int16(decVal);

filename = 'Audio_12b.wav';
audiowrite(filename,decVal_song,Fs,'BitsPerSample',16,...
'Comment','This is my new audio file.');
figure();
spectrogram(decVal);
title('Espectrograma de 12Bits');

%% Creación de filtros pasa bajas
%filtro 1500khz
F_filtro_max = Fs/2;
Fcorte_deseada = 15000;
Fcorte = Fcorte_deseada/F_filtro_max;

f = [0 Fcorte Fcorte 1];
m = [1 1 0 0]; o = 50;
LPF_15k = fir2(o,f,m);
fvtool(LPF_15k);  % Espectro en frecuencia del filtro.

%y = conv(x,b); % filtrado de x con filtro b o usando:
%y = filter(b,1,x]); % filtrado de x con filtro b, %función filter
%y = fftfilt(b,x);   % A more efficient FIR filtering for large operands

%%
%filtro 4kHz
Fcorte_deseada = 4000;
Fcorte = Fcorte_deseada/F_filtro_max;
f = [0 Fcorte Fcorte 1];
LPF_4k = fir2(o,f,m);
fvtool(LPF_4k);  % Espectro en frecuencia del filtro.


%%
%filtro 1kHz
Fcorte_deseada = 1000;
Fcorte = Fcorte_deseada/F_filtro_max;
f = [0 Fcorte Fcorte 1];
o = 70;
LPF_1k = fir2(o,f,m);
fvtool(LPF_1k);  % Espectro en frecuencia del filtro.
%%
%CREACION DE LAS CANCIONES CON FILTRO

Filtered_song = conv(y,LPF_15k); % filtrado de x con filtro b o usando:
filename = 'Audio_filtrado_15k.wav';
audiowrite(filename,Filtered_song,Fs,'BitsPerSample',16,...
'Comment','This is my new audio file.');
[X,Fs] = audioread('Audio_filtrado_15k.wav');
figure();
spectrogram(X);
title('Espectrograma de canción filtrada 15k');

Filtered_song = conv(y,LPF_4k); % filtrado de x con filtro b o usando:
filename = 'Audio_filtrado_4k.wav';
audiowrite(filename,Filtered_song,Fs,'BitsPerSample',16,...
'Comment','This is my new audio file.');
[X,Fs] = audioread('Audio_filtrado_4k.wav');
figure();
spectrogram(X);
title('Espectrograma de canción filtrada 4k');


Filtered_song = conv(y,LPF_1k); % filtrado de x con filtro b o usando:
filename = 'Audio_filtrado_1k.wav';
audiowrite(filename,Filtered_song,Fs,'BitsPerSample',16,...
'Comment','This is my new audio file.');
[X,Fs] = audioread('Audio_filtrado_1k.wav');
figure();
spectrogram(X);
title('Espectrograma de canción filtrada 1k');



