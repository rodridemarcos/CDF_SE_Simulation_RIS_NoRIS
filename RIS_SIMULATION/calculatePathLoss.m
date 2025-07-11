function PL = calculatePathLoss(scenario, frequency, distance, d_bp, hBS, hMS, isLOS)
    % Función auxiliar para calcular la pérdida de propagación según el escenario y LoS
    
    fc_GHz = frequency / 1e9; % Convertir la frecuencia a GHz
    d3D = distance; % Suponiendo que la distancia 3D es igual a la 2D por simplicidad aquí
    
    if strcmp(scenario, 'UMa')
        if isLOS == 1
            if distance >= 10 && distance <= d_bp
                PL = 28 + 22*log10(d3D) + 20*log10(fc_GHz);
            elseif distance > d_bp && distance <= 200 % Distancia para UMa
                PL = 28 + 40*log10(d3D) + 20*log10(fc_GHz) - 9*log10((d_bp)^2 + (hBS - hMS)^2);
            else
                PL = 100; % Valor alto para indicar fuera de rango
            end
        else % NLOS
            PL = 13.54 + 39.08*log10(d3D) + 20*log10(fc_GHz) - 0.6*(hMS - 1.5);
        end
    elseif strcmp(scenario, 'UMi')
        if isLOS == 1
            if distance >= 10 && distance <= d_bp
                PL = 32.4 + 21*log10(d3D) + 20*log10(fc_GHz);
            elseif distance > d_bp && distance <= 200 % Distancia para UMi
                PL = 32.4 + 40*log10(d3D) + 20*log10(fc_GHz) - 9.5*log10((d_bp)^2 + (hBS - hMS)^2);
            else
                PL = 100; % Valor alto para indicar fuera de rango
            end
        else % NLOS
            PL = 35.3*log10(d3D) + 22.4 + 21.3*log10(fc_GHz) - 0.3*(hMS - 1.5);
        end
    else
        error('Escenario no válido. Debe ser UMa o UMi.');
    end
end
