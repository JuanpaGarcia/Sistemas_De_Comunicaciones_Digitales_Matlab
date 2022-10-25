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
Rs = ((2*B)/Rs)-1;
Tp = 1/Rs; 
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
 energy = Tp;
 type = 'srrc';
 mp = round(Tp/Ts); %mp
 SRRC = rcpulse(beta,D,Tp,Ts,type,energy);

lena512 = imread('lena.tif');
pixels = 128; 
lenarec=lena512((284-(pixels-1)):284, (350-(pixels-1)):350); 
imshow(uint8(lenarec)) 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector
V_bit_pixels = b(1:pixels*pixels*8); %8 because 8bit pixel

V_bit_pixels_polar = zeros(1,numel(V_bit_pixels)*mp);

counter = 0;
for i= 0 : numel(V_bit_pixels)-1
    if V_bit_pixels(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_bit_pixels_polar(counter*i+1) = value;
    counter = mp;
end

srrc_PNRZ_sig = conv(V_bit_pixels_polar, SRRC);

if plots == 1
    bits_2_plot = 16;
    plot(srrc_PNRZ_sig(51:51+mp*bits_2_plot));
end 


start = floor(numel(SRRC)/2);
sampled_srrc_PNRZ_sig = srrc_PNRZ_sig(start:mp:end);

if plots == 1
    scatterplot(sampled_srrc_PNRZ_sig);
    title('Bits plot');
    %analize PSD
    figure();
    pwelch(srrc_PNRZ_sig,500,300,500,Fs,'power');
    title('PSD bit stream');
    %eyediagram 
    eyediagram(srrc_PNRZ_sig, 3*mp);
end 

%%
%Exercise 5
%Reception

%channel 
%Low pass filter generation, communications channel

f=[0 0.6 0.6 1];
m=[1 1 0 0];
ford=60;
filter_1 = fir2(ford,f,m);

%signal transmited
signal_filtered = conv(filter_1,srrc_PNRZ_sig);

%match filter
match_f = fliplr(SRRC);

recover_signal = conv(match_f,signal_filtered);

start_recovery = ford/2 + start*2;

recovery_signal = recover_signal(start_recovery:mp:end);
treshhold = 0;

if plots == 1
    figure();
    pwelch(recover_signal,500,300,500,Fs,'power');
    title('PSD bit stream');
    %eyediagram 
    eyediagram(recover_signal(start_recovery:mp*(128*128/4)), 3*mp);
    scatterplot(recovery_signal);
end 


PNRZ_recovery_bits_rec = zeros(1,numel(V_bit_pixels));
PNRZ_recovery_bits_rec((recovery_signal > treshhold)) = 1; 

%calculate BER

bits_error_PNRZ_hs = xor(V_bit_pixels, PNRZ_recovery_bits_rec(1:numel(V_bit_pixels)));
error_PNRZ = sum(bits_error_PNRZ_hs);
Bit_error_rate_PNRZ = (error_PNRZ/numel(V_bit_pixels)) * 100

%%
%Recover the image

%Half sine
recuperado = zeros(pixels,pixels,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PNRZ_recovery_bits_rec(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end

imshow(uint8(recuperado)) %recovered image
title('Image recovered');

%%

%%
%Exercise 6
%Transmition 

 beta = 0.5;
 Fs = 96000;
 Ts = 1/Fs;
 Rs = 9600;
 D = 10;
 Rb = Rs;
 Tp = 1/Rb; 
 energy = Tp;
 type = 'srrc';
 mp = round(Tp/Ts); %mp
 SRRC = rcpulse(beta,D,Tp,Ts,type,energy);

lena512 = imread('lena.tif');
pixels = 128; 
lenarec=lena512((284-(pixels-1)):284, (350-(pixels-1)):350); 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector
V_bit_pixels = b(1:pixels*pixels*8); %8 because 8bit pixel

V_bit_pixels_polar = zeros(1,numel(V_bit_pixels)*mp);

counter = 0;
for i= 0 : numel(V_bit_pixels)-1
    if V_bit_pixels(i+1) == 0
        value = -1;
    else
        value = 1;
    end
    V_bit_pixels_polar(counter*i+1) = value;
    counter = mp;
end

srrc_PNRZ_sig = conv(V_bit_pixels_polar, SRRC);

if plots == 3
    bits_2_plot = 16;
    plot(srrc_PNRZ_sig(51:51+mp*bits_2_plot));
end 


start = floor(numel(SRRC)/2);
sampled_srrc_PNRZ_sig = srrc_PNRZ_sig(start:mp:end);

if plots == 1
    scatterplot(sampled_srrc_PNRZ_sig);
    title('Bits plot');
    %analize PSD
    figure();
    pwelch(srrc_PNRZ_sig,500,300,500,Fs,'power');
    title('PSD bit stream');
    %eyediagram 
    eyediagram(srrc_PNRZ_sig(start:mp*2000), 3*mp);
end 
%%
%Reception

%channel 
%Low pass filter generation, communications channel

f=[0 0.6 0.6 1];
m=[1 1 0 0];
ford=60;
filter_1 = fir2(ford,f,m);

%signal transmited
signal_filtered = conv(filter_1,srrc_PNRZ_sig);

%match filter
match_f = fliplr(SRRC);

recover_signal = conv(match_f,signal_filtered);

start_recovery = ford/2 + start*2;

recovery_signal = recover_signal(start_recovery:mp:end);
treshhold = 0;

if plots == 1
    figure();
    pwelch(recover_signal,500,300,500,Fs,'power');
    title('PSD bit stream');
    %eyediagram 
    eyediagram(recover_signal(start_recovery:mp*(128*128/4)), 3*mp);
    scatterplot(recovery_signal);
end 


PNRZ_recovery_bits_rec = zeros(1,numel(V_bit_pixels));
PNRZ_recovery_bits_rec((recovery_signal > treshhold)) = 1; 

%calculate BER

bits_error_PNRZ_hs = xor(V_bit_pixels, PNRZ_recovery_bits_rec(1:numel(V_bit_pixels)));
error_PNRZ = sum(bits_error_PNRZ_hs);
Bit_error_rate_PNRZ = (error_PNRZ/numel(V_bit_pixels)) * 100


%Recover the image
%%
%Half sine
recuperado = zeros(pixels,pixels,'uint32'); %allocate memory
%load and convert values into a matrix
counter = 8; %counter variable
for i = 1 : pixels
    for j = 1: pixels
        recuperado(j,i) = bi2de(PNRZ_recovery_bits_rec(counter-7:counter),'right-msb');
        counter = counter +8;
    end
end

imshow(uint8(recuperado)) %recovered image
title('Image recovered');















