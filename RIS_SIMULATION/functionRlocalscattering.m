function R = functionRlocalscattering(M,angle,ASDdeg,antennaSpacing)
% Genera la matriz de correlación espacial para un canal con dispersión local.
%
% Entradas:
%   M             : Número de antenas.
%   angle         : Ángulo de llegada de la señal (radianes).
%   ASDdeg        : Desviación estándar angular en grados.
%   antennaSpacing: Espaciado entre antenas (número de longitudes de onda).
%
% Salidas:
%   R             : Matriz de correlación espacial (M x M).

ASD = ASDdeg * pi / 180; % Convertir a radianes
diff_angles = repmat(0:(M-1), M, 1) - repmat((0:(M-1))', 1, M);

R = exp(-((pi * antennaSpacing * diff_angles * sin(angle)).^2) / (2 * ASD^2));

% Rotar la matriz de correlación según el ángulo nominal de salida
varphi = (0:M-1)' * 2 * pi * antennaSpacing * sin(angle);
rotationVector = exp(1i * varphi);
R = rotationVector * rotationVector' .* R;

end
