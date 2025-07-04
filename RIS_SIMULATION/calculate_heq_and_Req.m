function [H_eq, R_eq] = calculate_heq_and_Req(HMean_BS_UE, H_BS_UE, H_BS_RIS, H_RIS_UE, theta_all, nbrOfRealizations, M_BS, K, L)

% Calcula el canal equivalente H_eq y su matriz de covarianza R_eq
% para el enlace descendente (BS -> RIS -> UE y directo BS -> UE)

% Entradas:
%   H_BS_UE         : Canal directo BS -> UE (M_BS x nbrOfRealizations x K)
%   H_BS_RIS        : Canal BS -> RIS (N_ris x M_BS x nbrOfRealizations)
%   H_RIS_UE        : Canal RIS -> UE (K x nbrOfRealizations x N_ris)
%   theta_all       : Matriz de fases de la RIS (N_ris x N_ris x nbrOfRealizations)
%   nbrOfRealizations : Número de realizaciones de canal
%   M_BS            : Número de antenas en la BS
%   K               : Número de usuarios
%   L               : Número de caminos (o celdas si aplica)
%   N_ris           : Número de elementos de la RIS

% Salidas:
%   H_eq            : Canal equivalente (M_BS x nbrOfRealizations x K)
%   R_eq            : Matriz de covarianza (M_BS x M_BS x K)
%   HMean_eq        : Media del canal equivalente (M_BS x K)

% Inicialización de matrices
H_eq_aux = zeros(K, M_BS, nbrOfRealizations, L, L);
R_eq = zeros(M_BS, M_BS, K, L, L);
%HMean_eq = zeros(M_BS, K);

for t = 1:nbrOfRealizations
    for j = 1:L
        for l = 1:L
            % Canal reflejado: UE <- RIS <- BS
            % h_reflected: (K x N_ris) * (N_ris x N_ris) * (N_ris x M_BS) = (K x M_BS)
            h_reflect = squeeze(H_RIS_UE(:, t, :, j, l))' * ...
                        theta_all(:, :, t, l) * ...
                        squeeze(H_BS_RIS(:, :, t, j ,l)).';
            
            % Canal equivalente: directo + reflejado
            % H_BS_UE(:,t,:) = (M_BS x K) → se transpone para suma con h_reflect (K x M_BS)'
            H_eq_aux(:, :, t, j, l) = squeeze(H_BS_UE(:, t, :, j, l)).' + h_reflect;
        
    

            % Reordenar dimensiones: (K x M_BS x realizaciones) → (M_BS x realizaciones x K)
            H_eq = permute(H_eq_aux, [2, 3, 1, 4, 5]);
            %H_eq = H_eq;
        
            % Cálculo de R_eq (matriz de covarianza) y HMean_eq (media)
            for k = 1:K
                % Extraer todas las realizaciones para el usuario k (M_BS x realizaciones)
                H_k = squeeze(H_eq(:, t, k, j, l));
                
                % Media del canal
                %HMean_eq(:, k) = mean(H_k, 2);
                
                % Calcular covarianza: R = E{h * h^H}
                %for l = 1:size(H_k, 2)
                    R_eq(:, :, k, j, l) = R_eq(:, :, k, j, l) + (H_k - HMean_BS_UE(:, k,j,l)) * (H_k - HMean_BS_UE(:, k,j,l))';
                %end
            end
        end
    end
end
% Promediar sobre el número de realizaciones
R_eq = R_eq / nbrOfRealizations;

end