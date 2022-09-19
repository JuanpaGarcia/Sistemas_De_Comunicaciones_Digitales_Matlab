%Miranda Elizabeth Dávila Velarde

%Tarea 4

%Ejercicio 1



%%

%Ejercicio 2 y 3

Fs=96000;
Ts=1/Fs;
mp=20;

%Unipolar NRZ
pbase = rectwin(mp); 
wvtool(pbase);
s = zeros(1,numel(bits)*mp);
s(1:mp:end) = bits;   % Impulse Train
figure;
stem(s(1:mp*16));
title('Impulse Train for Unipolar NRZ Stem');
figure;
plot(s(1:mp*16));
title('Impulse Train for Unipolar NRZ Plot');
xUNRZ = conv(pbase,s); % Pulse Train
figure;
plot(xUNRZ(1:mp*16));
title('Convolución Unipolar NRZ');
figure;
pwelch(xUNRZ,[],[],[],Fs,'power');  % PSD of Unipolar NRZ

%%
%Polar NRZ
pbase = rectwin(mp);
wvtool(pbase);
s1 = bits;
s1(s1==0) = -1; 
s = zeros(1,numel(s1)*mp);
s(1:mp:end) = s1;
figure;
stem(s(1:mp*16));
title('Impulse Train for Polar NRZ Stem');
figure;
plot(s(1:mp*16));
title('Impulse Train for Polar NRZ Plot');
xPNRZ = conv(pbase,s);
figure;
plot(xPNRZ(1:mp*16));
title('Convolución Polar NRZ');
figure;
pwelch(xPNRZ,500,300,500,Fs,'power');

%%
%Polar RZ
pbase2 = rectwin(mp/2);
wvtool(pbase2);
s2=bits;
s2(s2==0) = -1; 
s = zeros(1,numel(s2)*mp);
s(1:mp:end) = s2;
figure;
stem(s(1:mp*16));
title('Impulse Train for Polar RZ Stem');
figure;
plot(s(1:mp*16));
title('Impulse Train for Polar RZ Plot');
xPRZ = conv(pbase2,s);
figure;
plot(xPRZ(1:mp*16));
title('Convolución Polar RZ');
figure;
pwelch(xPRZ,500,300,500,Fs,'power');


%%
%Bipolar NRZ
pbase = rectwin(mp);
wvtool(pbase);
s2=bits;
lastbit=1;
for i=1:length(bits)
  if bits(i)==1
    s2(i)=lastbit;
    lastbit=-lastbit;
  end
end
s = zeros(1,numel(s2)*mp);
s(1:mp:end) = s2;
figure;
stem(s(1:mp*16));
title('Impulse Train for Bipolar NRZ Stem');
figure;
plot(s(1:mp*16));
title('Impulse Train for Bipolar NRZ Plot');
xBNRZ = conv(pbase2,s);
figure;
plot(xBNRZ(1:mp*16));
title('Convolución Bipolar NRZ');
figure;
pwelch(xBNRZ,500,300,500,Fs,'power');

%%
%Manchester
pbase2 = rectwin(mp);
wvtool(pbase2);
s3=bits;
offset = 1;

for i = 0:length(bits)-1
    if bits(i+offset) == 1
        s3(i+offset:(0.5+i+offset)) = 1;
        s3((i+0.5+offset):i+1+offset) = -1;
    else
        s3((i+offset):(0.5+i+offset)) = -1;
        s3((i+0.5+offset):(i+1+offset)) = 1;
    end
end

figure;
stem(s3(1:mp*16));
title('Impulse Train for Manchester Stem');
figure;
plot(s3(1:mp*16));
title('Impulse Train for Manchester Plot');
xM = conv(pbase2,s3);
figure;
plot(xM(1:mp*16));
title('Convolución Manchester');
figure;
pwelch(xM,500,300,500,Fs,'power');

