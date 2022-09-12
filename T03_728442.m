clear all;
close all;
clc 

%Generar pulsos en matlab
t = [0 0.1 0.2 0.3 0.4 0.5]; 
x = zeros(1,numel(t)); 
x(find( (t>=0.2) & (t<=0.4) )) = 1; 

%%
%Ejercicio 1
fs = 100;
t1 = 0:(1/fs):1-(1/fs);
x1 = zeros(1,numel(t1));
x1(find( (t1>=0.4) & (t1<=0.6) )) = 1; 
%mirar wtool
wvtool(x1)
%%
%Segundo pulso 
x2 = zeros(1,numel(t1));
x2(find( (t1>=0.0) & (t1<=0.2) )) = 1;
wvtool(x2)
%%
%Segundo ejercicio
x3 = zeros(1,numel(t1));
x3(find( ((t1>=0.2) & (t1<=0.4)) | ((t1>=0.7) & (t1<=0.9) ) )) = 1;
wvtool(x3)

%%
%tercer ejercicio
x4 = zeros(1,numel(t1));
x4(find( (t1>=0.2) & (t1<=0.4) )) = 1;
x4(find( (t1>=0.7) & (t1<=0.9) )) = -1;
wvtool(x4)

%%
%4to ejercicio
x5 = zeros(1,numel(t1));
x6 = zeros(1,numel(t1));
x7 = zeros(1,numel(t1));
x8 = zeros(1,numel(t1));

x5(find( (t1>=0.45) & (t1<=0.55) )) = 1;
x6(find( (t1>=0.475) & (t1<=0.525) )) = 1;
x7(find( (t1>=0.3) & (t1<=0.7) )) = 1;
x8(find( (t1>=0.2) & (t1<=0.8) )) = 1;

wvtool(x5)






