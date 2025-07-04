function PL = calculatePathLoss_UMi(frequency, distance, d_bp, hBS, hMS, isLOS)
    % Helper function to calculate path loss based on UMi scenario and LoS
    
    fc_GHz = frequency / 1e9; % Convert frequency to GHz
    d3D = distance; % Assuming 3D distance is the same as 2D for simplicity here
    
    if isLOS == 1
        if distance >= 10 && distance <= d_bp
            PL = 32.4 + 21*log10(d3D) + 20*log10(fc_GHz);
        elseif distance > d_bp && distance <= 500 % Distancia para UMi (ajustado a 500m)
            PL = 32.4 + 40*log10(d3D) + 20*log10(fc_GHz) - 9.5*log10((d_bp)^2 + (hBS - hMS)^2);
        else
            PL = 100; % Valor alto para indicar fuera de rango
        end
    else % NLOS
        PL = 35.3*log10(d3D) + 22.4 + 21.3*log10(fc_GHz) - 0.3*(hMS - 1.5);
    end
end