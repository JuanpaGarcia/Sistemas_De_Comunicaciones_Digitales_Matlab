clear all;
close all;
clc 

%Tarea 2    

[y,Fs] = audioread('spring_HIFI.wav');

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



