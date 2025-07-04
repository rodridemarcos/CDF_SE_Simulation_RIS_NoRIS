function [R_BS_UE, HMean_BS_UE, channelGaindB_BS_UE, ricianFactor_BS_UE, probLOS_BS_UE, worstUserIndex_per_cell, bestUsersIndices_per_cell] = functionExampleSetup(L,K,M,ASDdeg,scenario,frequency, LoS, seed)
% Genera la configuración de ejemplo para la simulación, incluyendo la
% ubicación de las BSs y los UEs, y calcula las propiedades del canal
% (pérdidas de trayectoria, probabilidad de LoS, factor de Rice, correlación
% espacial y media del canal) para el enlace BS-UE.
% Este script está configurado para el escenario UMi con una distancia de 500 metros y sin considerar RIS.
%
% Entradas:
%   L           : Número de estaciones base (BSs).
%   K           : Número de usuarios (UEs) por BS.
%   M           : Número de antenas en cada BS.
%   ASDdeg      : Desviación estándar angular (Angular Standard Deviation) en grados.
%   scenario    : 'UMa' para Urbano Macro o 'UMi' para Urbano Micro. (NOTA: Este script solo soporta UMi)
%   frequency   : Frecuencia de operación en hercios.
%   LoS         : Indicador para la condición de línea de vista:
%                 0: NLOS forzado para todos los enlaces.
%                 1: LOS forzado para todos los enlaces.
%                 2: Probabilidad de LOS dependiente de la distancia.
%                 3: Caso especial definido por el usuario (9 LOS, 1 NLOS).
%   seed        : Semilla para el generador de números aleatorios (si > 0).
%
% Salidas:
%   R_BS_UE             : Matriz de correlación espacial de las antenas BS-UE (M x M x K x L x L).
%   HMean_BS_UE         : Vector de la media del canal BS-UE (M x K x L x L).
%   channelGaindB_BS_UE : Ganancia del canal BS-UE en dB (K x L x L).
%   ricianFactor_BS_UE  : Factor de Rice BS-UE (K x L x L).
%   probLOS_BS_UE       : Probabilidad de línea de vista BS-UE (K x L x L).
%   worstUserIndex_per_cell: Índice del usuario NLOS para cada celda (1 x L).
%   bestUsersIndices_per_cell: Celdas con índices de usuarios LOS para cada celda (1 x L cell array).

% Constantes
c = 3e8; % Velocidad de la luz

% Alturas de la BS y el MS (solo UMi)
hBS = 10; % metros (altura de la estación base en UMi)
hMS = 1.5; % metros (altura del terminal móvil)

% Calcular la distancia de punto de ruptura (breakpoint distance)
hBS_prime = hBS - 1; % Altura efectiva de la BS
hMS_prime = hMS - 1; % Altura efectiva del MS
d_bp = (4 * hBS_prime * hMS_prime * frequency) / c; % Distancia de punto de ruptura

% Parámetros del modelo de ubicación
squareLength = 500; % metros (lado del área de cobertura para UMi)
nbrBSsPerDim = sqrt(L); % Número de BSs por dimensión
minDistance = 35; % metros (distancia mínima permitida entre un UE y su BS)
antennaSpacing = 1/2; % Espaciado de antenas (en número de longitudes de onda)
interBSDistance = squareLength/nbrBSsPerDim;
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
% Para manejar el envoltorio de las celdas (wrap-around)
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

% Inicialización de matrices para BS-UE
UEpositions = zeros(K,L);
perBS = zeros(L,1);
R_BS_UE = zeros(M,M,K,L,L);
HMean_BS_UE = zeros(M,K,L,L);
channelGaindB_BS_UE = zeros(K,L,L);
probLOS_BS_UE = zeros(K,L,L);
ricianFactor_BS_UE = zeros(K,L,L);
maxDistLOS = 300; % metros (distancia máxima para considerar LOS en el modelo probabilístico)

% Salidas para índices de usuarios
worstUserIndex_per_cell = zeros(1, L);
bestUsersIndices_per_cell = cell(1, L);

if seed > 0
    rng(seed); % Inicializar el generador de números aleatorios
end

for l = 1:L % Iterar sobre la celda propia (donde el UE está sirviendo)
    % Asignar UEs a la celda
    while perBS(l)<K
        UEremaining = K-perBS(l);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;
        posXY = posXY(abs(posXY)>=minDistance); % Asegurar distancia mínima
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
    end

    % --- Calcular ganancias para BS-UE y encontrar peor/mejor usuario en la celda l ---
    channelGaindB_local_BS_UE = zeros(K,1);
    for k = 1:K
        distance_BS_UE_local = abs(UEpositions(k,l) - BSpositions(l));
        % Calcular PL para el enlace directo (asumimos LOS para este cálculo de ordenación)
        PL = calculatePathLoss_UMi(frequency, distance_BS_UE_local, d_bp, hBS, hMS, 1); % Siempre LOS para ordenación
        channelGaindB_local_BS_UE(k) = -PL;
    end
    [~, sortedIndices] = sort(channelGaindB_local_BS_UE);
    worstUserIndex_per_cell(l) = sortedIndices(1); % El usuario con la peor ganancia en la celda l
    bestUsersIndices_per_cell{l} = sortedIndices(2:K); % Los K-1 usuarios con mejor ganancia

    % --- Bucle para las BS interferentes (j) ---
    for j = 1:L % Iterar sobre la BS desde la que se recibe (propia o interferente)
        % --- Cálculos para el enlace BS-UE ---
        [distancesBSj,~] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        for k = 1:K
            % Probabilidad de LOS
            if LoS == 1
                probLOS_BS_UE(k,l,j) = 1;
            elseif LoS == 0
                probLOS_BS_UE(k,l,j) = 0;
            elseif LoS == 2 % Probabilidad de LOS dependiente de la distancia
                probLOS_BS_UE(k,l,j) = (rand(1) < ((maxDistLOS - distancesBSj(k)) ./ maxDistLOS));
                if probLOS_BS_UE(k,l,j) < 0
                    probLOS_BS_UE(k,l,j) = 0;
                end
                if distancesBSj(k) > maxDistLOS
                    probLOS_BS_UE(k,l,j) = 0;
                end
            elseif LoS == 3 % Caso especial: 9 LOS, 1 NLOS
                if j == l % Para la BS de la propia celda (l)
                    if k == worstUserIndex_per_cell(l)
                        probLOS_BS_UE(k, l, j) = 0; % El peor usuario es NLOS
                    else
                        probLOS_BS_UE(k, l, j) = 1; % El resto es LOS
                    end
                else % Para las BSs de las otras celdas (interferentes)
                    probLOS_BS_UE(k, l, j) = 0; % Todos los usuarios son NLOS
                end
            end
            
            % Calcular el Ricean Factor
            if probLOS_BS_UE(k,l,j)==1
                riceFactor_dB_BS_UE = 13 - 0.03 * distancesBSj(k);
                ricianFactor_BS_UE(k,l,j) = db2pow(riceFactor_dB_BS_UE);
            else
                ricianFactor_BS_UE(k,l,j) = 0; % Factor de Rice es 0 para NLOS
            end

            % Calcular las pérdidas de trayectoria BS-UE
            PL_BS_UE = calculatePathLoss_UMi(frequency, distancesBSj(k), d_bp, hBS, hMS, probLOS_BS_UE(k,l,j));
            channelGaindB_BS_UE(k,l,j) = -PL_BS_UE;
        end
        
        % --- Aplicar Shadowing (para cada enlace) ---
        for k = 1:K
            if probLOS_BS_UE(k,l,j)==1
                shadowing_BS_UE = 4*randn(1); % sigma_sf_LOS=4 dB
            else
                shadowing_BS_UE = 10*randn(1); % sigma_sf_NLOS=10 dB
            end
            channelGaindB_BS_UE(k,l,j) = channelGaindB_BS_UE(k,l,j) + shadowing_BS_UE;
        end

        % --- Calcular la media del canal y la matriz de correlación espacial ---
        for k=1:K
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,min(find(abs(UEpositions(k,l) - BSpositionsWrapped(j,:)) == distancesBSj(k)))));
            HMean_BS_UE(:,k,l,j)=(exp(1i*2*pi.*(0:(M-1))'*sin(angleBSj)*antennaSpacing));
            R_BS_UE(:,:,k,l,j) = functionRlocalscattering(M,angleBSj,ASDdeg,antennaSpacing);
        end
    end
end
end
