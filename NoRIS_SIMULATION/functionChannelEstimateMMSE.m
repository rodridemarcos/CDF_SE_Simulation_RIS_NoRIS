function [Hhat_MMSE, C_MMSE] = functionChannelEstimateMMSE(R,HMeanx,H,nbrOfRealizations,M,K,L,p,f,tau_p)
% Generación de estimaciones de canal MMSE
%
%This Matlab function was developed to generate simulation results to:
%
%Ozgecan Ozdogan, Emil Bjornson, Erik G. Larsson, “Massive MIMO with
%Spatially Correlated Rician Fading Channels,” IEEE Transactions on
%Communications, To appear.
%
%Download article: https://arxiv.org/abs/1805.07972
%
%This is version 1.0 (Last edited: 2019-02-01)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


% Remodelar los vectores promedio para obtener el mismo vector promedio para todas las realizaciones
HMean=reshape(repmat(HMeanx,nbrOfRealizations,1),M,nbrOfRealizations,K,L,L);

% Generar el patrón de pilotos
if f == 1
    pilotPattern = ones(L,1);
elseif f == 2 % Solo funciona en el ejemplo ejecutándose con sus 16 BSs
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 2; 1; 2; 1]);
elseif f == 4 % Solo funciona en el ejemplo ejecutándose con sus 16 BSs
    pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]);
elseif f == 16 % Solo funciona en el ejemplo ejecutándose con sus 16 BSs
    pilotPattern = (1:L)';
elseif f == K % Agregado para manejar el caso f=K para el patrón de pilotos
    pilotPattern = (1:L)';
end

% Almacenar matriz identidad de tamaño M x M
eyeM = eye(M);

% Generar realizaciones de ruido normalizado
Np = sqrt(0.5)*(randn(M,nbrOfRealizations,K,L,f) + 1i*randn(M,nbrOfRealizations,K,L,f));

% Preparar almacenamiento para las estimaciones MMSE del canal
Hhat_MMSE = zeros(M,nbrOfRealizations,K,L,L);

% Preparar almacenamiento para las matrices de covarianza de error MMSE
C_MMSE = zeros(M, M, K, L, L); % Inicializar C_MMSE

% Recorrer todas las celdas
for j = 1:L

    % Recorrer todos los grupos de pilotos f
    for g = 1:f

        % Extraer las celdas que pertenecen al grupo de pilotos g
        groupMembers = find(g==pilotPattern)';

        % Calcular la señal de piloto procesada para todos los UEs que usan estos pilotos, según (5)
        yp = sqrt(p)*tau_p*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(tau_p)*Np(:,:,:,j,g);
        yMean=sqrt(p)*tau_p*sum(HMean(:,:,:,g==pilotPattern,j),4);
        % Recorrer todos los UEs
        for k = 1:K

            % Calcular la matriz que se invierte en el estimador MMSE
            PsiInv = (p*tau_p*sum(R(:,:,k,groupMembers,j),4) + eyeM);

            % Recorrer las celdas en el grupo de pilotos g
            for l = groupMembers

                % Calcular la estimación MMSE del canal entre la BS l y el UE k en
                % la celda j usando (6) en el Lema 1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = HMean(:,:,k,l,j) + sqrt(p)*RPsi*(yp(:,:,k)-yMean(:,:,k));

                % Calcular la matriz de covarianza de error MMSE usando (8) en el Lema 1
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);

            end

        end

    end

end
end
