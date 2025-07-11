function R = functionRlocalscattering(M,angleofdeparture,ASDdeg,antennaSpacing)
% Calcula la matriz de correlación espacial basada en el modelo de dispersión local.
%
% ENTRADAS:
% M                 = Número de antenas
% angleofdeparture  = Ángulo nominal de salida (en radianes)
% ASDdeg            = Desviación estándar angular (en grados)
% antennaSpacing    = Espaciado de antenas (en número de longitudes de onda)
%
% SALIDA:
% R                 = Matriz de correlación espacial M x M
%
%
%This Matlab function was developed to generate simulation results to:
%
%Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017),
%"Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency",
%Foundations and Trends in Signal Processing: Vol. 11, No. 3-4,
%pp. 154-655. DOI: 10.1561/2000000093.
%
%For further information, visit: https://www.massivemimobook.com
%
%This is version 1.0 (Last edited: 2017-11-04)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


% Convertir ASD de grados a radianes
ASD = ASDdeg*pi/180;

% Preparar matriz para la salida
R = zeros(M,M);

% Recorrer todos los pares de antenas
for m1 = 1:M
    for m2 = 1:M

        % Calcular la distancia entre las antenas
        delta_m = (m1-m2)*antennaSpacing;

        % Calcular la correlación espacial de acuerdo con (2.15) del monográfico
        R(m1,m2) = besselj(0,2*pi*delta_m*sin(ASD/2));

    end
end

% Rotar la matriz de correlación de acuerdo al ángulo nominal de salida
varphi = (0:M-1)'*2*pi*antennaSpacing*sin(angleofdeparture);
rotationVector = exp(1i*varphi);
R = rotationVector*rotationVector'.*R;
