function [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_UL(Hhat,C,R,tau_c,tau_p,nbrOfRealizations,M,K,L,p, N0)
%Compute UL SE for different receive combining schemes using Theorem 4.1.
%
%INPUT:
%Hhat              = M x nbrOfRealizations x K x L x L matrix with the MMSE
%                      channel estimates
%C                 = M x M x K x L x L matrix with estimation error
%                      correlation matrices when using MMSE estimation
%R                 = M x M x K x L x L matrix with spatial correlation
%                      matrices
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per BS
%K                 = Number of UEs per cell
%L                 = Number of BSs and cells
%p                 = Uplink transmit power per UE (same for everyone)
%N0                = Noise variance (in Watts)
%
%OUTPUT:
%SE_MR      = K x L matrix where element (k,l) is the uplink SE of UE k in
%             cell l achieved with MR combining
%SE_RZF     = Same as SE_MR but with RZF combining
%SE_MMMSE   = Same as SE_MR but with M-MMSE combining
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


%Store identity matrices of different sizes
eyeK = eye(K);
eyeM = eye(M);

%Compute the pre-log factor (normalized with number of channel realizations)
%assuming only uplink transmission
prelogFactor = (tau_c-tau_p)/(tau_c*nbrOfRealizations);

%Prepare to store simulation results
SE_MR = zeros(K,L);

if nargout > 1
    SE_RZF = zeros(K,L);
end

if nargout > 2
    SE_MMMSE = zeros(K,L);
    
    %Compute sum of all estimation error correlation matrices at every BS
    C_totM = reshape(p*sum(sum(C,3),4),[M M L]);
    
end

%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Go through all cells
    for j = 1:L
        
        %Extract channel estimate realizations from all UEs to BS j
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M K*L]);
        
        %Compute MR combining in (4.11)
        V_MR = Hhatallj(:,K*(j-1)+1:K*j);
        
        if nargout > 1 %Compute RZF combining in (4.9)
            V_RZF = (p * V_MR) / (p * (V_MR' * V_MR) + eyeK);
        end
        
        if nargout > 2 %Compute M-MMSE combining in (4.7)
            V_MMMSE = (p * (Hhatallj * Hhatallj') + C_totM(:,:,j) + N0*eyeM) \ (p * V_MR);
        end
                
        
        %Go through all UEs in cell j
        for k = 1:K
            
            %%MR combining
            v = V_MR(:,k); %Extract combining vector
            
            %Compute numerator and denominator of instantaneous SINR in (4.3)
            numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
            denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+N0*eyeM)*v - numerator;
            
            %Compute instantaneous SE for one channel realization
            SE_MR(k,j) = SE_MR(k,j) + prelogFactor*real(log2(1+numerator/denominator));
            
                    
            
            %%RZF combining
            if nargout > 1
                
                v = V_RZF(:,k); %Extract combining vector
                
                %Compute numerator and denominator of instantaneous SINR in (4.3)
                numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
                denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+N0*eyeM)*v - numerator;
                
                %Compute instantaneous SE for one channel realization
                SE_RZF(k,j) = SE_RZF(k,j) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
            
            %%M-MMSE combining
            if nargout > 2
                
                v = V_MMMSE(:,k); %Extract combining vector
                
                %Compute numerator and denominator of instantaneous SINR in (4.3)
                numerator = p*abs(v'*Hhat(:,n,k,j,j))^2;
                denominator = p*sum(abs(v'*Hhatallj).^2) + v'*(C_totM(:,:,j)+N0*eyeM)*v - numerator;
                
                %Compute instantaneous SE for one channel realization
                SE_MMMSE(k,j) = SE_MMMSE(k,j) + prelogFactor*real(log2(1+numerator/denominator));
                
            end
            
        end
        
    end
    
end
end