clear all;close all;
% se al muy similar a la anal gica��
% Cosenoidal de frecuencia f Hz
f=10;
Ts = 1/1000;
t = 0:Ts:1;
amp = 1;
x = amp*cos(2*pi*t*f);
% se al con muestreo natural 19 veces m s lento (53Hz)��
tm = t(1:19:end);
xm = cos(2*pi*tm*f);
% se al muestreada con sample-and-hold�
tsh = 19;
xs = zeros(1,numel(t));
for i=1:numel(t)
if( rem(i,tsh)==1 )
tmp = x(i);
end
xs(i) = tmp;
end
% Proceso de Cuantificaci n�
M = 65536;
int = (max(xs)-min(xs))/M;
m = (min(xs)+int/2):int:(max(xs)-int/2);
xq = zeros(1,length(t));
for i=1:length(t)
[tmp k] = min(abs(xs(i)-m));
xq(i) = m(k);
end
% diferencia
xd = xs - xq;
% graficas
figure(1)
plot(t,x)  %analog signal
title('Imagen figura analogica vs sampling');
xlabel('periodo') 
ylabel('amplitud') 
hold on
stem(tm,xm,'r','filled') %natural sampling
figure(2)
plot(t,xs)  %sample and hold
title('Imagen Sample and hold');
xlabel('periodo') 
ylabel('amplitud')
figure(3)
plot(t,xq) %quantified
title('Imagen cuantificación');
xlabel('periodo') 
ylabel('amplitud')
figure(4)
plot(t,xd) %difference
title('Imagen diferencias ');
xlabel('periodo') 
ylabel('amplitud')

figure(5);
tiledlayout(2,2)

nexttile
plot(t,x)  %analog signal
title('Imagen figura analogica vs sampling');
xlabel('periodo') 
ylabel('amplitud') 
hold on
stem(tm,xm,'r','filled') %natural sampling

nexttile
plot(t,xs)  %sample and hold
title('Imagen Sample and hold');
xlabel('periodo') 
ylabel('amplitud')

nexttile
plot(t,xq) %quantified
title('Imagen cuantificación');
xlabel('periodo') 
ylabel('amplitud')

nexttile
plot(t,xd) %difference
title('Imagen diferencias ');
xlabel('periodo') 
ylabel('amplitud')

%Potencia de la se al: �
Px=x*x'/numel(x); %var(x)
ex = xs-xq; %Forma 1
%%*****************Forma 2
%  b=log2(M);
%  xq= round(x*(2^(b-1)))/(2^(b-1)); %Se al Cuantizada�
%  ex = x-xq; 
%*************************************
Pn=ex*ex'/numel(ex);
SQNR_xdB = 10*log10(Px/Pn) %Simulado
SQNR_Teorico2= 6.02*(log2(M)) + 1.76 % Para una se al sinusoidal