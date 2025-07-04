function [R_out,HMean_out,H,H_Rayleigh] = functionChannelGeneration( RNormalized_in,HMeanNormalized_in,channelGaindB,ricianFactor,probLOS,K_in,L_in,M_in,nbrOfRealizations)
% R_out: Matriz de covarianza de la componente NLOS del canal, escalada.
%        Sus dimensiones son size(RNormalized_in,1) x size(RNormalized_in,2) x K_in x L_in x L_in.
% HMean_out: Matriz de la componente LOS del canal (valor medio), escalada.
%            Sus dimensiones son M_in x K_in x L_in x L_in.
% H: Realizaciones del canal total (LOS + NLOS).
%    Sus dimensiones son M_in x nbrOfRealizations x K_in x L_in x L_in.
% H_Rayleigh: Realizaciones del canal de Rayleigh (solo componente NLOS).
%             Sus dimensiones son M_in x nbrOfRealizations x K_in x L_in x L_in.

% RNormalized_in: Matriz de covarianza espacial normalizada (entrada).
% HMeanNormalized_in: Matriz del valor medio del canal normalizado (entrada).
% channelGaindB: Ganancia del canal en dB.
% ricianFactor: Factor K de Rician.
% probLOS: Probabilidad de Line-of-Sight (LOS).
% K_in: Número de usuarios (o elementos RIS que actúan como usuarios en un enlace).
% L_in: Número de células o dimensión de sectorización/interferencia.
% M_in: Número de antenas del array (primera dimensión del canal H).
% nbrOfRealizations: Número de realizaciones de Monte Carlo.

% Preparar almacenamiento para ganancias de canal por componente
channelGain_LOS=zeros(K_in,L_in,L_in);    % Ganancia para el componente LOS
channelGain_NLOS=zeros(K_in,L_in,L_in);  % Ganancia para el componente NLOS

% Inicializar la matriz del valor medio del canal (componente LOS)
HMean_out=zeros(M_in,K_in,L_in,L_in);

% Obtener las dimensiones de la matriz de covarianza de entrada para inicializar R_out
size_R_dim1 = size(RNormalized_in, 1);
size_R_dim2 = size(RNormalized_in, 2);

% Inicializar la matriz de covarianza escalada (componente NLOS)
R_out=zeros(size_R_dim1,size_R_dim2,K_in,L_in,L_in); 

% Iterar sobre usuarios/RIS elementos (k), células (l) y (j) para escalar matrices de covarianza y vectores medios
for k=1:K_in    % This 'k' is now L (4)
    for l=1:L_in  % This 'l' is now N_ris (64)
        for j=1:L_in  % This 'j' is also N_ris (64)

            % Accessing probLOS, ricianFactor, channelGaindB:
            % They are assumed to be (K_in, L_in, 1) or (K_in, L_in)
            if probLOS(k,l,j)==1 % <-- This now means probLOS(L_idx, N_ris_idx, 1)
                channelGain_LOS(k,l,j)= sqrt(ricianFactor(k,l,j)/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
                channelGain_NLOS(k,l,j)=sqrt(1/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
            else
                channelGain_LOS(k,l,j)= 0;
                channelGain_NLOS(k,l,j)=db2pow(channelGaindB(k,l,j));
            end

            % HMean_out and R_out assignments:
            HMean_out(:,k,l,j)=sqrt(channelGain_LOS(k,l,j))*HMeanNormalized_in(:,k,l,j);
            R_out(:,:,k,l,j)=channelGain_NLOS(k,l,j)*RNormalized_in(:,:,k,l,j);

        end
    end
end

% Generar las realizaciones del canal
% Generar variables aleatorias no correlacionadas para la componente Rayleigh
W = (randn(M_in,nbrOfRealizations,K_in,L_in,L_in)+1i*randn(M_in,nbrOfRealizations,K_in,L_in,L_in));

% Inicializar las matrices para las realizaciones del canal
H=zeros(M_in,nbrOfRealizations,K_in,L_in,L_in);      % Canal total
H_Rayleigh=zeros(M_in,nbrOfRealizations,K_in,L_in,L_in); % Canal de Rayleigh (solo NLOS)

% Expandir HMean_out para que tenga una dimensión para cada realización
% HMean_out: M_in x K_in x L_in x L_in
% HMeanx_expanded: M_in x nbrOfRealizations x K_in x L_in x L_in
HMeanx_expanded = permute(repmat(HMean_out, [1 1 1 1 nbrOfRealizations]), [1 5 2 3 4]);


% Iterar sobre usuarios/RIS elementos (k), células (j), y células (l) para generar realizaciones
for j = 1:L_in
    for l = 1:L_in
        for k = 1:K_in
            % Calcular la raíz cuadrada de la matriz de covarianza para la componente Rayleigh
            Rsqrt = sqrtm(R_out(:,:,k,j,l));
            
            % Generar la realización del canal total (LOS + NLOS)
            H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l) + HMeanx_expanded(:,:,k,j,l);
            
            % Generar la realización del canal de Rayleigh (solo NLOS)
            H_Rayleigh(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l);
        end
    end
end

end