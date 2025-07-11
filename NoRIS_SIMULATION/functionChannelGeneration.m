function [R,HMean,H,H_Rayleigh] = functionChannelGeneration( RNormalized,HMeanNormalized,channelGaindB,ricianFactor,probLOS,K,L,M,nbrOfRealizations)
% Escala las matrices de covarianza normalizadas y los vectores promedio por el
% factor de Rician y la ganancia de canal.
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

% Preparar variables para almacenar
channelGain_LOS=zeros(K,L,L);
channelGain_NLOS=zeros(K,L,L);
HMean=zeros(M,K,L,L);
R=zeros(M,M,K,L,L);

% Recorrer todos los UEs y aplicar las ganancias de canal a las matrices de
% correlación espacial y a los vectores promedio
for k=1:K
    for l=1:L
        for j=1:L
            
            if probLOS(k,l,j)==1 % Existe camino LoS, Factor de Rician ~= 0
                channelGain_LOS(k,l,j)= sqrt(ricianFactor(k,l,j)/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
                channelGain_NLOS(k,l,j)=sqrt(1/(ricianFactor(k,l,j) +1 ))*db2pow(channelGaindB(k,l,j));
            else  % Caso puramente NLoS
                channelGain_LOS(k,l,j)= 0;
                channelGain_NLOS(k,l,j)=db2pow(channelGaindB(k,l,j));
            end
            % Operación de escalado
            HMean(:,k,l,j)=sqrt(channelGain_LOS(k,l,j))*HMeanNormalized(:,k,l,j);
            R(:,:,k,l,j)=channelGain_NLOS(k,l,j)*RNormalized(:,:,k,l,j);
            
        end
    end
end


% Generar las realizaciones de canal
% Generar variables aleatorias no correlacionadas
W = (randn(M,nbrOfRealizations,K,L,L)+1i*randn(M,nbrOfRealizations,K,L,L));

% Preparar variables para almacenar las realizaciones de canal
H=zeros(M,nbrOfRealizations,K,L,L);
H_Rayleigh=zeros(M,nbrOfRealizations,K,L,L);

% Remodelar los vectores promedio para obtener el mismo vector promedio para todas las realizaciones
HMeanx=reshape(repmat(HMean,nbrOfRealizations,1),M,nbrOfRealizations,K,L,L);

% Recorrer todos los UEs y generar las realizaciones de canal
for j = 1:L
    
    for l = 1:L
        
        for k = 1:K
            Rsqrt = sqrtm(R(:,:,k,j,l));
            H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l) + HMeanx(:,:,k,j,l);
            H_Rayleigh(:,:,k,j,l) = sqrt(0.5)*Rsqrt*W(:,:,k,j,l);
        end
        
    end
    
end

end

   
