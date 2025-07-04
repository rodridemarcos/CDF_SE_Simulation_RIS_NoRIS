function [Hhat_MMSE, C_MMSE] = functionChannelEstimateMMSE(R, HMeanx, H, nbrOfRealizations, M, K, L, p, f, tau_p)
                               
% Genera estimaciones MMSE del canal para un sistema Massive MIMO

% Expandimos el vector de medias para tener una media igual en todas las realizaciones
HMean = reshape(repmat(HMeanx, nbrOfRealizations, 1), M, nbrOfRealizations, K, L, L);

% Generamos el patrón de pilotos (asignación de grupo de piloto para cada celda)
if f == 1
    pilotPattern = ones(L,1); % Todos usan el mismo piloto
elseif f >= L
    pilotPattern = (1:L)';    % Cada estación base tiene un piloto diferente
else
    % Asignación cíclica de pilotos para adaptarse a cualquier L y f
    pilotPattern = mod((0:L-1), f) + 1;
end

% Matriz identidad M x M
eyeM = eye(M);

% Generamos ruido normalizado complejo (distribución gaussiana circular)
Np = sqrt(0.5)*(randn(M, nbrOfRealizations, K, L, f) + 1i*randn(M, nbrOfRealizations, K, L, f));

% Inicializamos almacenamiento para las estimaciones MMSE
Hhat_MMSE = zeros(M, nbrOfRealizations, K, L, L);

% Inicializamos almacenamiento para las matrices de covarianza de error MMSE
C_MMSE = zeros(M, M, K, L, L);

% Recorremos todas las celdas receptoras (estaciones base)
for j = 1:L

    % Recorremos todos los grupos de pilotos
    for g = 1:f

        % Celdas que usan el mismo grupo de piloto g
        groupMembers = find(g == pilotPattern)';

        % Calculamos la señal piloto procesada según la ecuación (5)
        yp = sqrt(p)*tau_p*sum(H(:,:,:,g == pilotPattern, j), 4) + sqrt(tau_p)*Np(:,:,:,j,g);
        yMean = sqrt(p)*tau_p*sum(HMean(:,:,:,g == pilotPattern, j), 4);

        % Recorremos todos los usuarios
        for k = 1:K

            % Matriz que será invertida en el estimador MMSE
            PsiInv = (p*tau_p*sum(R(:,:,k,groupMembers,j), 4) + eyeM);

            % Recorremos todas las celdas en el grupo de piloto g
            for l = groupMembers

                % Estimación MMSE del canal según (6) del Lema 1
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = HMean(:,:,k,l,j) + sqrt(p)*RPsi*(yp(:,:,k) - yMean(:,:,k));

                % Matriz de covarianza de error MMSE según (8) del Lema 1
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p*tau_p*RPsi*R(:,:,k,l,j);

            end
        end
    end
end
end
