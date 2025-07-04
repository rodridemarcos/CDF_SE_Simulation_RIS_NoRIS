function Theta = calculate_theta(H_BS_UE, H_BS_RIS, H_RIS_UE, N_ris, nbrOfRealizations, L, p, worstuser)

% Calcula las fases óptimas de la RIS en downlink
%
% Entradas:
%   H_BS_UE        : Canal directo BS -> UE (M x realizaciones x K x L x L)
%   H_BS_RIS       : Canal BS -> RIS (N_ris x realizaciones x M x L x L)
%   H_RIS_UE       : Canal RIS -> UE (N_ris x realizaciones x K x L x L)
%   N_ris          : Nº de elementos RIS
%   nbrOfRealizations : Nº de realizaciones de canal
%   K              : Nº de UEs por celda
%   L              : Nº de celdas
%   p              : SNR en escala lineal
%   fase           : Es la theta o fase de cada metaátomo de la RIS
%   Hs'            : Es la conjugada traspuesta de Hs
%   Hn'            : Es la conjugada traspuesta de Hn
%
% Salida:
%   Theta          : Matriz de fases óptimas de la RIS (N_ris x N_ris x realizaciones x L)

Theta = zeros(N_ris, N_ris, nbrOfRealizations, L);

for t = 1:nbrOfRealizations
    for l_cell = 1:L

        k = worstuser(l_cell);  % Usuario seleccionado para optimizar la RIS en esta celda

        % Canal directo BS → UE (M x 1)
        Hs = squeeze(H_BS_UE(:, t, k, l_cell, l_cell));  % M x 1

        % Canal BS → RIS (Ht): entrada a la RIS (M x N_ris)
        Hrx = squeeze(H_BS_RIS(:, :, t, l_cell, l_cell));  % M x N_ris

        % Canal RIS → UE (Hr): salida desde la RIS (N_ris x 1)
        Htx = squeeze(H_RIS_UE(:, t, k, l_cell, l_cell));  % N_ris x 1

        M = size(Htx, 2);  % Nº de antenas en la BS

        % Inicializar vector de fase aleatoria para RIS
        fase = exp(1i * 2 * pi * rand(N_ris, 1));

        for n = 1:N_ris

            % Canal efectivo con todas las reflexiones menos la del elemento n
            % Hr' = transpuesta conjugada de Hr (1 x N_ris)
            % diag(psiVec) * Ht → N_ris x M
            % Resultado Hn: M x 1
            Hn = Hs + Hrx * diag(fase) * Htx - fase(n) * Hrx(:,n) * Htx(n,:);

            % Vector bn (M x 1): depende del canal efectivo y del canal RIS→UE y BS→RIS en el elemento n
            bn = p * Hn* Htx(n,:)';  % M x 1

            % Matriz An (M x M)
            % Hn * Hn' → M x M (producto de vector columna por su transpuesta conjugada)
            % Hr(n) * Ht(n,:) → 1 x M, su transpuesta conjugada: M x 1
            An = eye(M) + p * (Hn * Hn') + p * (Hrx(:,n) * Htx(n,:)) * (Hrx(:,n) * Htx(n,:))';
            
            warning('off', 'all');%esta línea elimina los warnings del calculate theta

            % Actualizar fase del elemento n usando ángulo del producto escalar complejo
            fase(n) = exp(-1i * angle(bn' * (An \ (Hrx(:,n)))));

        end

        % Guardar matriz diagonal de fases óptimas
        Theta(:,:,t,l_cell) = diag(fase);

    end
end
    warning('on', 'all'); %reactiva los warnings que puedan saltar
end