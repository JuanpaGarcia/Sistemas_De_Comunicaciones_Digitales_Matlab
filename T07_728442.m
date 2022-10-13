clear all;
close all;
clc 

%Tarea 7

%%
mp = 10;
Fs = 96000;

lena512 = imread('lena.tif');
imshow(uint8(lena512)) 
lenarec=lena512(243:284,309:350); 
imshow(uint8(lenarec)) 

b=de2bi(lenarec,8); 
b=b'; 
bits=b(:);   % Bits vector

V_bit = b(1:42*42*8);

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
n = 0:mp-1;
w = pi/mp;
hs = sin(w*n);
pr = ones(1,mp);

% Polar NRZ LineCode rectangular
Polar_NRZ_sig_rec = conv(pr ,V_bit_polar);

% Polar NRZ LineCode halfsine
Polar_NRZ_sig_hs = conv(hs ,V_bit_polar);

%change power to 1
Polar_NRZ_sig_rec = sqrt(1/((sum(Polar_NRZ_sig_rec.^2))/numel(Polar_NRZ_sig_rec))).*Polar_NRZ_sig_rec;

%change power to 1
Polar_NRZ_sig_hs = sqrt(1/((sum(Polar_NRZ_sig_hs.^2))/numel(Polar_NRZ_sig_hs))).*Polar_NRZ_sig_hs;
%%
plot(Polar_NRZ_sig_rec);
title('Señal PNRZ rectangular');

figure();
plot(Polar_NRZ_sig_hs);
title('Señal PNRZ half sine');

figure();
eyediagram(Polar_NRZ_sig_rec,2*mp);
title('Eyediagrma rectangular');

figure();
eyediagram(Polar_NRZ_sig_hs,2*mp);
title('Eyediagrma half sine');

%%
%Ejercicio 3
%Low pass filter generation, communications channel
%F 0.3042
f=[0 0.4 0.4 1];
m=[1 1 0 0];
ford=60;
filter_1 = fir2(ford,f,m);
%%
%transmit the signal in the communication channel
Signal_filtered_Polar_NRZ_rec = conv(Polar_NRZ_sig_rec, filter_1);
plot(Signal_filtered_Polar_NRZ_rec)
title('Filtered signal Polar NRZ rectangular');


Signal_filtered_Polar_NRZ_hs = conv(Polar_NRZ_sig_hs, filter_1);
figure();
plot(Signal_filtered_Polar_NRZ_hs)
title('Filtered signal Polar NRZ half sine');

delay = ford/2;

figure();
eyediagram(Signal_filtered_Polar_NRZ_rec,2*mp);
title('Eyediagrma rectangular');

figure();
eyediagram(Signal_filtered_Polar_NRZ_hs,2*mp);
title('Eyediagrma half sine');

%%
%obtener ancho de banda 
% fs/mp
pwelch(Polar_NRZ_sig_hs);
title('Line coding bandwidth');

%%
%change power to 1
Signal_filtered_Polar_NRZ_hs = sqrt(1/((sum(Signal_filtered_Polar_NRZ_hs.^2))/numel(Signal_filtered_Polar_NRZ_hs))).*Signal_filtered_Polar_NRZ_hs;
plot(Signal_filtered_Polar_NRZ_hs)
title('Normalized Filtered signal Polar NRZ half sine');

%check power pow = (sum(Unipolar_NRZ_Sig.^2))/numel(Unipolar_NRZ_Sig);

%change power to 1
Signal_filtered_Polar_NRZ_rec = sqrt(1/((sum(Signal_filtered_Polar_NRZ_rec.^2))/numel(Signal_filtered_Polar_NRZ_rec))).*Signal_filtered_Polar_NRZ_rec;
figure();
plot(Signal_filtered_Polar_NRZ_rec)
title('Normalized Filtered signal Polar NRZ rectangular');

%%
%Ejercicio 5
%apply match filter

match_filtered_PNRZ_rec = conv(Signal_filtered_Polar_NRZ_rec, fliplr(pr));


match_filtered_PNRZ_hs = conv(Signal_filtered_Polar_NRZ_hs, fliplr(hs));

eyediagram(match_filtered_PNRZ_rec,2*mp);
title('Eyediagrma rectangular');

figure();
eyediagram(match_filtered_PNRZ_hs,2*mp);
title('Eyediagrma half sine');
























