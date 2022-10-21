% rcpulse.m
%
% Autor:  Luis Miguel Bazdresch Sierra
%         Departamento de Electrónica, Sistemas en Informática
%         Universidad ITESO
%         http://iteso.mx/~miguelbaz
%
% Cálculo de pulsos coseno alzado y raíz de coseno alzado
%
% argumentos de entrada:
%
%    beta = factor de rolloff
%    D = duración del pulso como número de intervalos de símbolo
%    Tp = intervalo de símbolo (cada cuánto tiempo se envía un nuevo símbolo)
%    Ts = intervalo de muestreo
%    type = tipo de pulso a generar (es un string)
%        'srrc' - square-root raised cosine
%        'rc'   - raised cosine (valor por default)
%    energy = energía del pulso, default = 1
%
% salidas:
%
%    p = las muestras del pulso
%    t = el vector de tiempo correspondiente al pulso
%
%
% referencias:
%  Telecommunication breakdown, páginas 217 y 225
%  http://www.dsplog.com/2008/04/22/raised-cosine-filter-for-transmit-pulse-shaping/
%  http://en.wikipedia.org/wiki/Raised-cosine_filter
%  http://en.wikipedia.org/wiki/Root-raised-cosine_filter

function [p t] = rcpulse(beta,D,Tp,Ts,type,energy)

if nargin < 6
	energy = 1;
end

if nargin == 4
	type = 'rc';
end

%% para evitar problemas numéricos con beta
if beta==0
	beta=1e-12;
end;

% checar que Tp sea múltiplo entero de Ts
if abs(round(Tp/Ts) - Tp/Ts) > 1e-6
	error('Error: Tp debe ser múltiplo de Ts');
end

% checar que D sea par
if mod(D,2) == 1
	error('Error: D debe ser par');
end

%%% vector de tiempo (no causal)
t = -D*Tp/2:Ts:D*Tp/2;

if strcmp(type, 'srrc')
	% square-root raised cosine (si beta==0, regresa un sinc)
	% se calcula en tres partes: tiempo negativo, 0, y positivo
	t1 = -D*Tp/2:Ts:-Ts;
	t2 = Ts:Ts:D*Tp/2;

	x1 = pi*(1-beta)*t1/Tp;
	x2 = pi*(1+beta)*t1/Tp;
	x3 = 4*beta*t1/Tp;
	num = sin(x1)+x3.*cos(x2);
	den = sqrt(Tp).*(pi*t1./Tp).*(1-x3.^2);
	p1 = num./den;
	p1(find(abs(den)<1e-9))=(beta/sqrt(2*Tp))*...
        ((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/(4*beta)));

	p2 = 1/sqrt(Tp)*(1-beta+4*beta/pi);

	x1 = pi*(1-beta)*t2/Tp;
	x2 = pi*(1+beta)*t2/Tp;
	x3 = 4*beta*t2/Tp;
	num = sin(x1)+x3.*cos(x2);
	den = sqrt(Tp).*(pi*t2./Tp).*(1-x3.^2);
	p3 = num./den;
	p3(find(abs(den)<1e-9))=(beta/sqrt(2*Tp))*...
		((1+2/pi)*sin(pi/(4*beta))+(1-2/pi)*cos(pi/(4*beta)));

	p = [p1 p2 p3];

elseif strcmp(type, 'rc')
	% raised cosine (si beta==0, regresa un sinc)
	% sólo hay que cuidar el caso donde el denominador se hace cero
	den = 1-(2*beta*t/Tp).^2;
	c = cos(pi*beta*t/Tp)./den;
	c(find(abs(den)<1e-9)) = pi/4;
	p = (1/Tp).*sinc(t/Tp).*c;
end

% hacer que el pulso sea causal
t = t + D*Tp/2;

% normalización
en = trapz(t,p.*p);
p = sqrt(energy/en)*p;
end
