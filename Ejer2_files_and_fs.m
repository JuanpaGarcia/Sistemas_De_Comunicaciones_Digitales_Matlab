clear all
load handel.mat % Recuperar la señal “y” en valores [-1 1]
k=8;                                 % Cantidad de bits para Cuantizar
%% Cuantizar a entero y expresión en binario
swing = (2^k-1)/2;                   % Señal simétrica
xq_int = round(y*swing+swing);       % Convertir a entero en [0 2^b-1]
xq_bin = de2bi(xq_int,k,'left-msb'); % Convertir entero a binario
% xq_bin = decimalToBinaryVector(xq_int,k,'MSBFirst'); %Opción 1
% xq_bin = decimalToBinaryVector(xq_int,k,'LSBFirst'); %Opción 2

xq = (xq_int-swing)/swing; % Proceso inverso usando xq_int

decVal= binaryVectorToDecimal(xq_bin,'MSBFirst');
figure();
plot(y);
figure();
plot(xq_bin);
[valor My] = max(y);

DtoA_Vec = (decVal-swing)/swing;
DtoA_Vec(My)

