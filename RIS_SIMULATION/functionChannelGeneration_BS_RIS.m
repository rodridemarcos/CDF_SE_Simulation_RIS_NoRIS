function [R_out_Tx, R_out_Rx, HMean_out, H, H_Rayleigh] = functionChannelGeneration_BS_RIS(RNormalized_in_Rx, RNormalized_in_Tx, HMeanNormalized_in, channelGaindB, ricianFactor, probLOS, N_ris, L, M, nbrOfRealizations)
% R_out_Rx: Matriz de covarianza de la componente NLOS del canal para el RECEPTOR (RIS), escalada.
%           Sus dimensiones son N_ris x N_ris x N_ris_loop x L x L.
% R_out_Tx: Matriz de covarianza de la componente NLOS del canal para el TRANSMISOR (BS), escalada.
%           Sus dimensiones son Mmax x Mmax x N_ris_loop x L x L.
% HMean_out: Matriz de la componente LOS del canal (valor medio), escalada.
%            Sus dimensiones son M x N_ris x L x L.
% H: Realizaciones del canal total (LOS + NLOS).
%    Sus dimensiones son M x N_ris x nbrOfRealizations x N_ris_loop x L x L.
% H_Rayleigh: Realizaciones del canal de Rayleigh (solo componente NLOS).
%             Sus dimensiones son M x N_ris x nbrOfRealizations x N_ris_loop x L x L.
% RNormalized_in_Rx: Matriz de covarianza espacial normalizada del RECEPTOR (entrada N_ris x N_ris x N_ris_loop x L x L).
% RNormalized_in_Tx: Matriz de covarianza espacial normalizada del TRANSMISOR (entrada Mmax x Mmax x N_ris_loop x L x L).
% HMeanNormalized_in: Matriz del valor medio del canal normalizado (entrada).
% channelGaindB: Ganancia del canal en dB.
% ricianFactor: Factor N_ris de Rician.
% probLOS: Probabilidad de Line-of-Sight (LOS).
% N_ris: Número de entidades en el lado del TRANSMISOR del canal H (Mmax para BS-RIS).
% L: Número de células o dimensión de sectorización/interferencia.
% M: Número de entidades en el lado del RECEPTOR del canal H (N_ris para BS-RIS).
% nbrOfRealizations: Número de realizaciones de Monte Carlo.

% Prepara el almacenamiento para las ganancias del canal por componente
channelGain_LOS = zeros(N_ris, L, L);     % Ganancia para el componente LOS
channelGain_NLOS = zeros(N_ris, L, L);   % Ganancia para el componente NLOS

% Inicializa la matriz del valor medio del canal (componente LOS)
HMean_out = zeros(M, N_ris, L, L);

% Obtener dimensiones de las matrices de covarianza de entrada
size_R_Rx_dim1 = size(RNormalized_in_Rx, 1); % Esperado N_ris
size_R_Rx_dim2 = size(RNormalized_in_Rx, 2); % Esperado N_ris
size_R_Tx_dim1 = size(RNormalized_in_Tx, 1); % Esperado Mmax
size_R_Tx_dim2 = size(RNormalized_in_Tx, 2); % Esperado Mmax

% Inicializa las matrices de covarianza escaladas (componente NLOS)
% R_out_Rx se indexa por N_ris,l,j, pero para la dimensión 3 se usa N_ris_loop (N_ris)
R_out_Rx = zeros(size_R_Rx_dim1, size_R_Rx_dim2,  L, L); % N_ris x N_ris x L x L
R_out_Tx = zeros(size_R_Tx_dim1, size_R_Tx_dim2, L, L); % Mmax x Mmax x L x L

% Iterar sobre el índice del elemento RIS (N_ris_idx), y las celdas (l, j)
% Para el enlace BS-RIS, N_ris_loop es N_ris, L es L.
for l = 1:L % Este 'l' es el primer índice de celda (1 a L)
    for j = 1:L % Este 'j' es el segundo índice de celda (1 a L)
        % Acceder a probLOS, ricianFactor, channelGaindB:
        % Nota: channelGaindB, ricianFactor, probLOS tienen dimensiones (N_ris x L x L) o (L x L)
        % Asumimos que channelGaindB(l,j), ricianFactor(l,j) y probLOS(N_ris_idx,l,j) son los correctos.
        if probLOS(l, j) == 1 % Asegurarse de que probLOS tenga la dimensión N_ris_idx
            channelGain_LOS(l, j) = sqrt(ricianFactor(l, j) / (ricianFactor(l, j) + 1)) * db2pow(channelGaindB(l, j));
            channelGain_NLOS(l, j) = sqrt(1 / (ricianFactor(l, j) + 1)) * db2pow(channelGaindB(l, j));
        else
            channelGain_LOS(l, j) = 0;
            channelGain_NLOS(l, j) = db2pow(channelGaindB(l, j));
        end
        
        % Asegurarse de que HMeanNormalized_in se indexa correctamente
        HMean_out(:, :, l, j) = sqrt(channelGain_LOS( l, j)) * HMeanNormalized_in(:, :, l, j); % HMeanNormalized_in: N_ris x Mmax x L x L
        
        % Escalar las matrices de correlación
        R_out_Rx(:, :, l, j) = channelGain_NLOS(l, j) * RNormalized_in_Rx(:, :, l, j);
        R_out_Tx(:, :, l, j) = channelGain_NLOS(l, j) * RNormalized_in_Tx(:, :, l, j);
    end
end

% Generar realizaciones del canal
% H y H_Rayleigh dimensiones: (M x N_ris x nbrOfRealizations x N_ris_loop x L x L)
H = zeros(M, N_ris, nbrOfRealizations, L, L);
H_Rayleigh = zeros(M, N_ris, nbrOfRealizations, L, L);

% Generar variables aleatorias no correlacionadas para el componente Rayleigh
% W_uncorrelated será M x N_ris para cada realización y cada N_ris_elem, j, l.
W_uncorrelated = (randn(M, N_ris, nbrOfRealizations, L, L) + ...
                  1i * randn(M, N_ris, nbrOfRealizations, L, L));

% Expandir HMean_out para la suma
HMeanx_expanded = repmat(HMean_out, [1 1 nbrOfRealizations 1 1]);

% Iterar sobre los elementos/usuarios (N_ris_elem), y las dos dimensiones de celda (j y l) para generar realizaciones
for j = 1:L
    for l = 1:L
        
        % Calcular las raíces cuadradas de las matrices de covarianza para el componente Rayleigh
        % R_out_Rx(:,:,N_ris_elem,j,l) es N_ris x N_ris (correlación Rx para este grupo de elementos RIS y celdas)
        Rsqrt_Rx = sqrtm(R_out_Rx(:, :, l, j));
        % R_out_Tx(:,:,N_ris_elem,j,l) es Mmax x Mmax (correlación Tx para este grupo de elementos RIS y celdas)
        Rsqrt_Tx = sqrtm(R_out_Tx(:, :, l, j));
        
        % Segmento de W_uncorrelated para el N_ris_elem,j,l actual sobre todas las realizaciones
        % W_slice_for_mult será M x N_ris x nbrOfRealizations (N_ris x Mmax x nbrOfRealizations)
        W_slice_for_mult = squeeze(W_uncorrelated(:, :, :, l, j));
        
        for real_idx = 1:nbrOfRealizations
            % H(:,:,real_idx,N_ris_elem,j,l) será M x N_ris (N_ris x Mmax)
            % HMeanx_expanded(:,:,real_idx,N_ris_elem,j,l) es M x N_ris (N_ris x Mmax)
            % W_slice_for_mult(:,:,real_idx) es M x N_ris (N_ris x Mmax)
            
            % Generar la realización del canal total (LOS + NLOS)
            % Dimensiones: (N_ris x N_ris) * (N_ris x Mmax) * (Mmax x Mmax) = (N_ris x Mmax)
            H_Rayleigh(:,:,real_idx, l, j) = sqrt(0.5) * Rsqrt_Tx * W_slice_for_mult(:,:,real_idx) * Rsqrt_Rx;
            H(:,:,real_idx,l, j) = H_Rayleigh(:,:,real_idx, l, j) + HMeanx_expanded(:,:, (real_idx-1)*L+l, j);
        end
        
    end
end
end