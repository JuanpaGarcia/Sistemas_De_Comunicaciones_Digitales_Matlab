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
load 'lena512.mat' 
imshow(uint8(lena512)) 
lenarec=lena512(252:284,318:350); 
imshow(uint8(lenarec)) 
%%
b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Vector de bits concatenado 


















