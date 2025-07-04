function [R_BS_UE, HMean_BS_UE, channelGaindB_BS_UE, ricianFactor_BS_UE, probLOS_BS_UE, ...
          R_UE_RIS, HMean_UE_RIS, channelGaindB_UE_RIS, ricianFactor_UE_RIS, probLOS_UE_RIS, ...
          R_BS_RIS_BSAnt, R_BS_RIS_RISel, HMean_BS_RIS, channelGaindB_BS_RIS, ricianFactor_BS_RIS, probLOS_BS_RIS, ...
          worstUserIndex_per_cell, bestUsersIndices_per_cell] = ...
    functionExampleSetupRIS(L,K,M,N_ris,ASDdeg,scenario,frequency, LoS, seed)
% Genera la configuración de ejemplo para la simulación, incluyendo la
% ubicación de las BSs, los UEs y las RISs, y calcula las propiedades del canal
% (pérdidas de trayectoria, probabilidad de LoS, factor de Rice, correlación
% espacial y media del canal) para los enlaces BS-UE, UE-RIS y BS-RIS.
%
% Entradas:ca
%   L             : Número de estaciones base (BSs).
%   K             : Número de usuarios (UEs) por BS.
%   M             : Número de antenas en cada BS.
%   ASDdeg        : Desviación estándar angular (Angular Standard Deviation) en grados.
%   scenario      : 'UMa' para Urbano Macro o 'UMi' para Urbano Micro.
%   frequency     : Frecuencia de operación en hercios.
%   LoS           : Indicador para la condición de línea de vista:
%                     0: NLOS forzado para todos los enlaces.
%                     1: LOS forzado para todos los enlaces.
%                     2: Probabilidad de LOS dependiente de la distancia.
%                     3: Caso especial definido por el usuario (9 LOS, 1 NLOS).
%   seed          : Semilla para el generador de números aleatorios (si > 0).
%
% Salidas:
%   R_BS_UE               : Matriz de correlación espacial de las antenas BS-UE (M x M x K x L x L).
%   HMean_BS_UE           : Vector de la media del canal BS-UE (M x K x L x L).
%   channelGaindB_BS_UE   : Ganancia del canal BS-UE en dB (K x L x L).
%   ricianFactor_BS_UE    : Factor de Rice BS-UE (K x L x L).
%   probLOS_BS_UE         : Probabilidad de línea de vista BS-UE (K x L x L).
%
%   R_UE_RIS              : Matriz de correlación espacial de los elementos UE-RIS (N_ris x N_ris x K x L x L).
%   HMean_UE_RIS          : Vector de la media del canal UE-RIS (N_ris x K x L x L).
%   channelGaindB_UE_RIS  : Ganancia del canal UE-RIS en dB (K x L x L).
%   ricianFactor_UE_RIS   : Factor de Rice UE-RIS (K x L x L).
%   probLOS_UE_RIS        : Probabilidad de línea de vista UE-RIS (K x L x L).
%
%   R_BS_RIS_BSAnt        : Matriz de correlación espacial de las antenas BS-RIS (M x M x L x L).
%   R_BS_RIS_RISel        : Matriz de correlación espacial de los elementos BS-RIS (N_ris x N_ris x L x L).
%   HMean_BS_RIS          : Vector de la media del canal BS-RIS (M x N_ris x L x L).
%   channelGaindB_BS_RIS  : Ganancia del canal BS-RIS en dB (L x L).
%   probLOS_BS_RIS        : Probabilidad de línea de vista BS-RIS (L x L).
%   ricianFactor_BS_RIS   : Factor de Rice BS-RIS (L x L).
%   worstUserIndex_per_cell: Índice del usuario NLOS para cada celda (1 x L).
%   bestUsersIndices_per_cell: Celdas con índices de usuarios LOS para cada celda (1 x L cell array).

% Constantes
c = 3e8; % Velocidad de la luz

% Alturas de la BS y el MS
if strcmp(scenario, 'UMa')
    hBS = 25; % metros
    hMS = 1.5; % metros
    squareLength = 500; % Área de cobertura para UMa
else % UMi
    hBS = 10; % metros
    hMS = 1.5; % metros
    squareLength = 500; % Área de cobertura para UMi
end

% Calcular la distancia de punto de ruptura
hBS_prime = hBS - 1;
hMS_prime = hMS - 1;
d_bp = (4 * hBS_prime * hMS_prime * frequency) / c;

% Parámetros del modelo de ubicación
nbrBSsPerDim = sqrt(L);
minDistance = 35; % Distancia mínima de UE a BS
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

% Número de elementos en la RIS (asumimos el mismo que antenas de la BS)
%N_ris = M; 

% Altura de la RIS 
hRIS = 10; % metros (altura intermedia entre BS y UE)

% Inicialización de matrices para BS-UE
UEpositions = zeros(K,L);
perBS = zeros(L,1);
R_BS_UE = zeros(M,M,K,L,L);
HMean_BS_UE = zeros(M,K,L,L);
channelGaindB_BS_UE = zeros(K,L,L);
probLOS_BS_UE = zeros(K,L,L);
ricianFactor_BS_UE = zeros(K,L,L);
maxDistLOS = 300; % Distancia máxima para considerar LOS en el modelo 3GPP

% --- Nuevas inicializaciones para RIS ---
% Inicialización de posiciones de la RIS (se calcularán dentro del bucle L)
RISpositions = zeros(1,L); 

% Inicialización de matrices para UE-RIS
R_UE_RIS = zeros(N_ris,N_ris,K,L,L);
HMean_UE_RIS = zeros(N_ris,K,L,L); 
channelGaindB_UE_RIS = zeros(K,L,L); % Ganancia UE-RIS (K UEs, L celdas, L RIS)
probLOS_UE_RIS = zeros(K,L,L);
ricianFactor_UE_RIS = zeros(K,L,L);

% Inicialización de matrices para BS-RIS
R_BS_RIS_BSAnt = zeros(M,M,L,L); % Correlación en las antenas de la BS
R_BS_RIS_RISel = zeros(N_ris,N_ris,L,L); % Correlación en los elementos de la RIS 
HMean_BS_RIS = zeros(M,N_ris,L,L); % Media del canal (M antenas BS, N_ris elementos RIS, L BS, L RIS)
channelGaindB_BS_RIS = zeros(L,L); % Ganancia BS-RIS (L BS, L RIS)
probLOS_BS_RIS = zeros(L,L);
ricianFactor_BS_RIS = zeros(L,L);

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

    % --- Calcular ganancias para BS-UE y encontrar peor usuario en la celda l ---
    channelGaindB_local_BS_UE = zeros(K,1);
    for k = 1:K
        distance_BS_UE_local = abs(UEpositions(k,l) - BSpositions(l));
        
        % Calcular PL para el enlace directo (solo para determinar el peor usuario)
        PL = calculatePathLoss(scenario, frequency, distance_BS_UE_local, d_bp, hBS, hMS, 1); 
        channelGaindB_local_BS_UE(k) = -PL;
    end
    [~, sortedIndices] = sort(channelGaindB_local_BS_UE);
    worstUserIndex_per_cell(l) = sortedIndices(1); % El usuario con la peor ganancia en la celda l
    bestUsersIndices_per_cell{l} = sortedIndices(2:K); % Los K-1 usuarios con mejor ganancia

    % --- Posicionar la RIS justo en medio de la BS y el peor usuario (NLOS) de la celda l --
    RISpositions(l) = (BSpositions(l) + UEpositions(worstUserIndex_per_cell(l), l)) / 2;
end 

% --- Una vez que todas las RISpositions están definidas, podemos crear RISpositionsWrapped ---
% Esto asegura que todas las RISpositions(l) están inicializadas antes de repmat.
%RISpositionsWrapped = repmat(RISpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);


% --- Bucle principal para calcular los parámetros de canal para cada enlace ---
for l = 1:L % Iterar sobre la celda propia (UE src)
    for j = 1:L % Iterar sobre la BS/RIS desde la que se recibe/interfiere
        %% --- Cálculos para el enlace BS-UE (Directo) ---
        % Las distancias de BS j a los UEs en la celda l.
        % Se busca la distancia mínima entre el UE(k,l) y cualquier BS(j) envuelta.
        [distancesBSj,~] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        % Probabilidad de LOS para BS-UE
        if LoS == 1
            probLOS_BS_UE(:,l,j) = ones(K,1);
        elseif LoS == 0
            probLOS_BS_UE(:,l,j) = zeros(K,1);
        elseif LoS == 2 % Probabilidad de LOS dependiente de la distancia
            probLOS_BS_UE(:,l,j) = (rand(K,1) < ((maxDistLOS - distancesBSj) ./ maxDistLOS));
            probLOS_BS_UE(probLOS_BS_UE < 0) = 0; % Asegurar que la probabilidad no sea negativa
            probLOS_BS_UE(distancesBSj > maxDistLOS) = 0; % Si la distancia es mayor que maxDistLOS, probLOS es 0
        elseif LoS == 3 % Caso especial: 9 LOS, 1 NLOS (el peor usuario)
            if j == l % Para la BS de la propia celda (l)
                probLOS_BS_UE(worstUserIndex_per_cell(l), l, j) = 0; % El peor usuario es NLOS
                probLOS_BS_UE(bestUsersIndices_per_cell{l}, l, j) = 1; % El resto es LOS
            else % Para las BSs de las otras celdas (interferentes)
                probLOS_BS_UE(:, l, j) = 0; % Todos los usuarios son NLOS
            end
        else
            error('Valor de LoS no válido. Debe ser 0, 1, 2 o 3.'); 
        end
        
        % Calcular el Ricean Factor en función de la distancia BS-UE
        riceFactor_dB_BS_UE = 13 - 0.03 * distancesBSj;
        ricianFactor_BS_UE(:,l,j) = db2pow(riceFactor_dB_BS_UE);
        ricianFactor_BS_UE(probLOS_BS_UE(:,l,j) == 0, l, j) = 0; % Asegurar que el factor de Rice es 0 para NLOS
        
        for k = 1:K
            distance_BS_UE = abs(UEpositions(k,l) - BSpositions(j)); % Distancia directa UE k (celda l) a BS j
            % Calcular las pérdidas de trayectoria BS-UE
            PL_BS_UE = calculatePathLoss(scenario, frequency, distance_BS_UE, d_bp, hBS, hMS, probLOS_BS_UE(k,l,j));
            channelGaindB_BS_UE(k,l,j) = -PL_BS_UE;
        end

        % --- Cálculo para el enlace UE-RIS (para UE k en celda l y RIS j) ---
        % Las distancias de UEs en celda l a la RIS j.
        distance_UE_RIS = abs(UEpositions(:,l) - RISpositions(j)); 
        
        for k = 1:K
            % Probabilidad de LOS UE-RIS: solo el peor usuario con su propia RIS
            % (asumiendo que la RIS j es la RIS de la celda l)
            if j == l && k == worstUserIndex_per_cell(l)
                probLOS_UE_RIS(k,l,j) = 1; % Peor usuario con RIS de su celda tiene LOS
            else
                probLOS_UE_RIS(k,l,j) = 0; % El resto es NLOS
            end
            
            % Ricean Factor UE-RIS
            riceFactor_dB_UE_RIS = 13 - 0.03 * distance_UE_RIS(k); % Usamos la misma fórmula
            ricianFactor_UE_RIS(k,l,j) = db2pow(riceFactor_dB_UE_RIS);
            ricianFactor_UE_RIS(probLOS_UE_RIS(k,l,j) == 0, l, j) = 0; % Factor de Rice es 0 para NLOS
            
            % Pérdidas de trayectoria UE-RIS (asumimos modelo similar a UE-BS, pero con hRIS como altura de la RIS)
            PL_UE_RIS = calculatePathLoss(scenario, frequency, distance_UE_RIS(k), d_bp, hBS, hRIS, probLOS_UE_RIS(k,l,j)); % Usar hRIS
            channelGaindB_UE_RIS(k,l,j) = -PL_UE_RIS;
        end

        % --- Cálculo para el enlace BS-RIS (para BS l y RIS j) ---
        % Distancia de la BS l a la RIS j
        distance_BS_RIS = abs(BSpositions(l) - RISpositions(j)); 
        
        % Probabilidad de LOS BS-RIS: solo si la BS y la RIS son de la misma celda
        if j == l
            probLOS_BS_RIS(l,j) = 1;
        else
            probLOS_BS_RIS(l,j) = 0;
        end
        
        % Ricean Factor BS-RIS
        
        if probLOS_BS_RIS(l,j)==0
            ricianFactor_BS_RIS(l,j) = 0;
        else
            riceFactor_dB_BS_RIS = 13 - 0.03 * distance_BS_RIS; 
            ricianFactor_BS_RIS(l,j) = db2pow(riceFactor_dB_BS_RIS);
        end
        
        % Pérdidas de trayectoria BS-RIS (modelo de propagación BS-BS, con hRIS para la RIS)
        PL_BS_RIS = calculatePathLoss(scenario, frequency, distance_BS_RIS, d_bp, hBS, hRIS, probLOS_BS_RIS(l,j)); % Usar hRIS
        channelGaindB_BS_RIS(l,j) = -PL_BS_RIS;
        
        % --- Aplicar Shadowing (para cada enlace) ---
        % Shadowing para BS-UE
        for k = 1:K
            if probLOS_BS_UE(k,l,j)==1
                shadowing_BS_UE = 4*randn(1); % sigma_sf_LOS=4 dB
            else
                shadowing_BS_UE = 10*randn(1); % sigma_sf_NLOS=10 dB
            end
            channelGaindB_BS_UE(k,l,j) = channelGaindB_BS_UE(k,l,j) + shadowing_BS_UE;
        end
        
        % Shadowing para UE-RIS
        for k = 1:K
            if probLOS_UE_RIS(k,l,j)==1
                shadowing_UE_RIS = 4*randn(1);
            else
                shadowing_UE_RIS = 10*randn(1);
            end
            channelGaindB_UE_RIS(k,l,j) = channelGaindB_UE_RIS(k,l,j) + shadowing_UE_RIS;
        end
        % Shadowing para BS-RIS
        if probLOS_BS_RIS(l,j)==1
            shadowing_BS_RIS = 4*randn(1);
        else
            shadowing_BS_RIS = 10*randn(1);
        end
        channelGaindB_BS_RIS(l,j) = channelGaindB_BS_RIS(l,j) + shadowing_BS_RIS;

        % --- Cálculos de HMean y R ---
        % HMean y R para BS-UE
        for k=1:K
            % Encontrar la BS envuelta más cercana para calcular el ángulo correctamente
            % Esto se hace para asegurar que el ángulo se calcula con respecto a la BS 'real' más cercana
            [~, minIdx] = min(abs(UEpositions(k,l) - BSpositionsWrapped(j,:)));
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,minIdx));
            
            HMean_BS_UE(:,k,l,j)=(exp(1i*2*pi.*(0:(M-1))'*sin(angleBSj)*antennaSpacing));
            R_BS_UE(:,:,k,l,j) = functionRlocalscattering(M,angleBSj,ASDdeg,antennaSpacing);
        end
        % HMean y R para UE-RIS
        for k=1:K
            angleUE_RIS = angle(RISpositions(j) - UEpositions(k,l)); % Ángulo de la RIS vista desde el UE
            HMean_UE_RIS(:,k,l,j)=(exp(1i*2*pi.*(0:(N_ris-1))'*sin(angleUE_RIS)*antennaSpacing)); % N_ris elementos
            R_UE_RIS(:,:,k,l,j) = functionRlocalscattering(N_ris,angleUE_RIS,ASDdeg,antennaSpacing); % N_ris x N_ris
        end
        
        % HMean y R para BS-RIS
        angleBS_RIS = angle(RISpositions(j) - BSpositions(l)); % Ángulo de la RIS vista desde la BS (Tx)
        
        % HMean para BS-RIS (M antenas de BS, N_ris elementos de RIS)
        arrayResponse_BS = exp(1i*2*pi.*(0:(M-1))'*sin(angleBS_RIS)*antennaSpacing);
        % Para la RIS, el ángulo de incidencia es el ángulo del rayo desde la BS.
        % Si angleBS_RIS es el ángulo de la RIS vista desde la BS, entonces
        % la RIS 'recibe' este rayo. El ángulo de la respuesta del array de la RIS
        % es el ángulo del rayo incidente.
        arrayResponse_RIS = exp(1i*2*pi.*(0:(N_ris-1))'*sin(angleBS_RIS)*antennaSpacing); % Ángulo de la RIS visto desde la BS
        HMean_BS_RIS(:,:,l,j) = arrayResponse_BS * arrayResponse_RIS.'; % Media del canal entre BS y RIS
        
        % R para BS-RIS (M antenas de BS, N_ris elementos de RIS)
        R_BS_RIS_BSAnt(:,:,l,j) = functionRlocalscattering(M,angleBS_RIS,ASDdeg,antennaSpacing); % Correlación en la BS
        R_BS_RIS_RISel(:,:,l,j) = functionRlocalscattering(N_ris,angleBS_RIS,ASDdeg,antennaSpacing); % Correlación en la RIS
        
    end 
end 

end 