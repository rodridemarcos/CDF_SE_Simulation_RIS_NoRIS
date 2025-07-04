function [R,HMean,channelGaindB,ricianFactor,probLOS] = functionExampleSetup(L,K,M,ASDdeg,scenario,frequency, LoS, seed)
% Genera la configuración de ejemplo para la simulación, incluyendo la
% ubicación de las BSs y los UEs, y calcula las propiedades del canal
% (pérdidas de trayectoria, probabilidad de LoS, factor de Rice, correlación
% espacial y media del canal).
%
% Entradas:
%   L         : Número de estaciones base (BSs).
%   K         : Número de usuarios (UEs) por BS.
%   M         : Número de antenas en cada BS.
%   ASDdeg    : Desviación estándar angular (Angular Standard Deviation) en grados.
%   scenario  : 'UMa' para Urbano Macro o 'UMi' para Urbano Micro.
%   frequency : Frecuencia de operación en hercios.
%   LoS       : Indicador para la condición de línea de vista:
%               0: NLOS forzado para todos los enlaces.
%               1: LOS forzado para todos los enlaces.
%               2: Probabilidad de LOS dependiente de la distancia.
%               3: Caso especial definido por el usuario.
%   seed      : Semilla para el generador de números aleatorios (si > 0).
%
% Salidas:
%   R             : Matriz de correlación espacial de las antenas de la BS (M x M x K x L x L).
%   HMean         : Vector de la media del canal (M x K x L x L).
%   channelGaindB : Ganancia del canal en dB (K x L x L).
%   ricianFactor  : Factor de Rice (K x L x L).
%   probLOS       : Probabilidad de línea de vista (K x L x L).
% Constantes
c = 3e8; % Velocidad de la luz
% Alturas de la BS y el MS
if strcmp(scenario, 'UMa')
    hBS = 25; % metros
    hMS = 1.5; % metros
else % UMi
    hBS = 10; % metros
    hMS = 1.5; % metros
end
% Calcular la distancia de punto de ruptura
hBS_prime = hBS - 1;
hMS_prime = hMS - 1;
d_bp = (4 * hBS_prime * hMS_prime * frequency) / c;
% Parámetros del modelo de ubicación
if strcmp(scenario, 'UMi')
    squareLength = 500; % Área de cobertura para UMi
else
    squareLength = 500; % Área de cobertura para UMa
end
nbrBSsPerDim = sqrt(L);
minDistance = 35;
antennaSpacing = 1/2;
interBSDistance = squareLength/nbrBSsPerDim;
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
UEpositions = zeros(K,L);
perBS = zeros(L,1);
R = zeros(M,M,K,L,L);
HMean=zeros(M,K,L,L);
channelGaindB=zeros(K,L,L);
probLOS=zeros(K,L,L);
ricianFactor=zeros(K,L,L);
maxDistLOS = 300; % Distancia máxima para considerar LOS (puede ser ajustada)
if seed > 0
    rng(seed); % Inicializar el generador de números aleatorios
end
for l = 1:L
    while perBS(l)<K
        UEremaining = K-perBS(l);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;
        posXY = posXY(abs(posXY)>=minDistance);
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
    end
    channelGaindB_local = zeros(K,1);
    for k = 1:K
        distance = abs(UEpositions(k,l) - BSpositions(l));
        if strcmp(scenario, 'UMa')
            if distance >= 10 && distance <= d_bp
                PL = 28 + 22*log10(distance) + 20*log10(frequency/1e9);
            elseif distance > d_bp && distance <= 200
                PL = 28 + 40*log10(distance) + 20*log10(frequency/1e9) - 9*log10((d_bp)^2 + (hBS - hMS)^2);
            else
                PL = 100;
            end
        elseif strcmp(scenario, 'UMi')
            if distance >= 10 && distance <= d_bp
                PL = 32.4 + 21*log10(distance) + 20*log10(frequency/1e9);
            elseif distance > d_bp && distance <= 200
                PL = 32.4 + 40*log10(distance) + 20*log10(frequency/1e9) - 9.5*log10((d_bp)^2 + (hBS - hMS)^2);
            else
                PL = 100;
            end
        else
            error('Escenario no válido. Debe ser UMa o UMi.');
        end
        channelGaindB_local(k) = -PL;
    end
    [~, sortedIndices] = sort(channelGaindB_local);
    worstUserIndex = sortedIndices(1);
    bestUsersIndices = sortedIndices(2:K);

    for j = 1:L
        [distancesBSj,~] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        % Calcular la probabilidad de LoS
        if LoS == 1
            probLOS(:,l,j) = ones(K,1); % Todos los usuarios tienen LoS
        elseif LoS == 0
            probLOS(:,l,j) = zeros(K,1); % Todos los usuarios tienen NLoS
        elseif LoS == 2
           probLOS(:,l,j) = (rand(K,1) < ((maxDistLOS - distancesBSj) ./ maxDistLOS));
            probLOS(probLOS < 0) = 0; % Asegurar que la probabilidad no sea negativa
            probLOS(distancesBSj > maxDistLOS) = 0; % Si la distancia es mayor que maxDistLOS, probLOS es 0
        elseif LoS == 3
            if j == l % Para la BS de la propia celda
                probLOS(worstUserIndex, l, j) = 0; % El usuario con peor ganancia tiene LoS = 0
                probLOS(bestUsersIndices, l, j) = 1; % Los K-1 usuarios con mejor ganancia tienen LoS = 1
            else % Para las BSs de las otras celdas
                probLOS(:, l, j) = 0; % Todos los usuarios tienen LoS = 0
            end
        else
            error('Valor de LoS no válido. Debe ser 0, 1, 2 o 3.');
        end
        % Calcular el Ricean Factor en función de la distanciaBSj (en dB)
        riceFactor_dB = 13 - 0.03 * distancesBSj;
        riceFactor_linear = db2pow(riceFactor_dB);
        ricianFactor(:,l,j) = riceFactor_linear;
        for k = 1:K
            distance = abs(UEpositions(k,l) - BSpositions(j)); % Distancia directa
            % Calcular las pérdidas de trayectoria
            if strcmp(scenario, 'UMa')
                if probLOS(k,l,j) == 1
                    if distance >= 10 && distance <= d_bp
                        PL = 28 + 22*log10(distance) + 20*log10(frequency/1e9);
                    elseif distance > d_bp && distance <= 200 % Distancia para UMa
                        PL = 28 + 40*log10(distance) + 20*log10(frequency/1e9) - 9*log10((d_bp)^2 + (hBS - hMS)^2);
                    else
                        PL = 100; % Valor alto para indicar fuera de rango
                         probLOS(k,l,j) = 0; % Si está fuera de rango LOS, se considera NLOS
                    end
                else % NLOS
                    PL = 13.54 + 39.08*log10(distance) + 20*log10(frequency/1e9) - 0.6*(hMS -1.5);
                end
            elseif strcmp(scenario, 'UMi')
                 if probLOS(k,l,j) == 1
                    if distance >= 10 && distance <= d_bp
                        PL = 32.4 + 21*log10(distance) + 20*log10(frequency/1e9);
                    elseif distance > d_bp && distance <= 200 % Distancia para UMi
                        PL = 32.4 + 40*log10(distance) + 20*log10(frequency/1e9) - 9.5*log10((d_bp)^2 + (hBS - hMS)^2);
                    else
                        PL = 100; % Valor alto para indicar fuera de rango
                        probLOS(k,l,j) = 0; % Si está fuera de rango LOS, se considera NLOS
                    end
                else % NLOS
                    PL = 35.3*log10(distance) + 22.4 + 21.3*log10(frequency/1e9) - 0.3*(hMS - 1.5);
                end
            else
                error('Escenario no válido. Debe ser UMa o UMi.');
            end
            channelGaindB(k,l,j) = -PL; % La ganancia es la negativa de la pérdida en dB
        end

    end
    for k=1:K
         if probLOS(k,l,l)==1
            shadowing = 4*randn(1,1,L); % sigma_sf_LOS=4 dB
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        else
            shadowing = 10*randn(1,1,L); % sigma_sf_NLOS=10 dB
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        end
        while channelGainShadowing(l) < max(channelGainShadowing)
            if probLOS(k,l,l)==1
                shadowing = 4*randn(1,1,L);
                channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            else
                shadowing = 10*randn(1,1,L);
                channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            end
        end
        channelGaindB(k,l,:) = channelGainShadowing;
    end
    for j=1:L
        [~,whichpos] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        for k=1:K
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,whichpos(k)));
            HMean(:,k,l,j)=(exp(1i*2*pi.*(0:(M-1))*sin(angleBSj)*antennaSpacing));
            R(:,:,k,l,j) = functionRlocalscattering(M,angleBSj,ASDdeg,antennaSpacing);
        end
    end
end