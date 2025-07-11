function [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p)
% Calcula la SE en el enlace ascendente para diferentes esquemas de combinación usando el Teorema 4.1.
%
% ENTRADAS:
% Hhat              = Matriz M x nbrOfRealizations x K x L x L con las estimaciones de canal MMSE
% C                 = Matriz M x M x K x L x L con las matrices de correlación del error de estimación usando MMSE
% R                 = Matriz M x M x K x L x L con matrices de correlación espacial
% tau_c             = Longitud del bloque de coherencia
% tau_p             = Longitud de las secuencias piloto
% nbrOfRealizations = Número de realizaciones del canal
% M                 = Número de antenas por estación base (BS)
% K                 = Número de usuarios por celda (UE)
% L                 = Número de estaciones base y celdas
% p                 = Potencia de transmisión en el enlace ascendente por usuario (igual para todos)
%
% SALIDAS:
% SE_MR    = Matriz K x L donde el elemento (k,l) es la SE UL del usuario k en la celda l obtenida con combinación MR
% SE_RZF   = Igual que SE_MR pero con combinación RZF
% SE_MMMSE = Igual que SE_MR pero con combinación M-MMSE
% SE_ZF    = Igual que SE_MR pero con combinación ZF
% SE_SMMSE = Igual que SE_MR pero con combinación S-MMSE
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
%This is version 1.01 (Last edited: 2020-05-15)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.


% Almacenar matrices identidad de diferentes tamaños
eyeK = eye(K);
eyeM = eye(M);

% Calcular el factor pre-log (normalizado con el número de realizaciones del canal)
% asumiendo solo transmisión en el enlace ascendente
prelogFactor = (tau_c-tau_p)/(tau_c*nbrOfRealizations);

% Preparar espacio para almacenar los resultados de la simulación
SE_MR = zeros(K,L);

if nargout > 1
    SE_RZF = zeros(K,L);
end

if nargout > 2
    SE_MMMSE = zeros(K,L);
    
    % Calcular la suma de todas las matrices de correlación del error de estimación en cada BS
    C_totM = reshape(p*sum(sum(C,3),4),[M M L]);
    
end

%% Recorrer todas las realizaciones del canal
for n = 1:nbrOfRealizations
    
    % Recorrer todas las celdas
    for j = 1:L
        
        % Extraer las realizaciones de la estimación de canal de todos los UEs hacia la BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        % Calcular la combinación MR en (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        if nargout > 1 % Calcular la combinación RZF en (4.9)
            V_RZF = (p * V_MR) / (p * (V_MR' * V_MR) + eyeK);
        end
        
        if nargout > 2 % Calcular la combinación M-MMSE en (4.7)
            V_MMMSE = (p * (Hhatallj * Hhatallj') + C_totM(:,:,j) + eyeM) \ (p * V_MR);
        end
               
        
        % Recorrer todos los UEs en la celda j
        for k = 1:K
            
            %% Combinación MR
            v = V_MR(:,k); % Extraer el vector de combinación
            
            % Calcular numerador y denominador del SINR instantáneo en (4.3)
            numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
            denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - numerator;
            
            % Calcular la SE instantánea para una realización del canal
            SE_MR(k,j) = SE_MR(k,j) + prelogFactor*real(log2(1+numerator/denominator));
            
                 
            
            %% Combinación RZF
            if nargout > 1
                
                v = V_RZF(:,k); % Extraer el vector de combinación
                
                % Calcular numerador y denominador del SINR instantáneo en (4.3)
                numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
                denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - numerator;
                
                % Calcular la SE instantánea para una realización del canal
                SE_RZF(k,j) = SE_RZF(k,j) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
            
            %% Combinación M-MMSE
            if nargout > 2
                
                v = V_MMMSE(:,k); % Extraer el vector de combinación
                
                % Calcular numerador y denominador del SINR instantáneo en (4.3)
                numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
                denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+eyeM)*v - numerator;
                
                % Calcular la SE instantánea para una realización del canal
                SE_MMMSE(k,j) = SE_MMMSE(k,j) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
        end
        
    end
    
end

