%%main_RIS_No_RIS
close all;
clear;
clc;
%% --- ENTRADAS CONFIGURABLES POR EL USUARIO --- %%
frequencies = [2e9,8e9,15e9,28e9];    % Frecuencias (Hz) para la simulación (fijo en 15 GHz como se solicitó)
LOS_settings = 3;      % Configuraciones de Línea de Vista (LOS) (9 LOS y 1 NLOS)
precoders = ["MR", "MMMSE"];      % Esquemas de precodificación a simular [MR, MMMSE, RZF]
ris_options = "RIS";        % Escenarios a simular: con RIS o sin RIS (fijo en RIS como se solicitó)

%% --- PARÁMETROS FIJOS DE SIMULACIÓN --- %%
L = 4;                 % Número de celdas (Estaciones Base - BS)
K = 10;                % Número de usuarios (Equipos de Usuario - UE) por celda
M = 40;                % Número de antenas en cada BS (M = Mmax)
Mmax = M;              % Alias para el número máximo de antenas en la BS
N_ris = 50;            % Número de elementos en el RIS [50,100,200]
nbrOfSetups = 100;       % Número de configuraciones de red para promediar
nbrOfRealizations = 100; % Número de realizaciones de canal por configuración (para promediar dentro de una configuración)
B = 20e6;              % Ancho de banda del sistema (Hz)
p = 0.1;               % Potencia de transmisión por UE (W)
noiseFigure = 7;       % Figura de ruido del receptor (dB)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure; % Varianza del ruido en dBm
noiseVarianceW = db2pow(noiseVariancedBm - 30); % Varianza del ruido en Watts
tau_c = 200;           % Longitud de coherencia de canal (símbolos)
tau_p = K;             % Longitud de la secuencia piloto (símbolos)
f = 1;                 % Factor de reutilización de secuencias piloto
ASDdeg = 10;           % Dispersión angular de salidas/llegadas (grados)
scenario = 'UMi';      % Escenario de propagación (por ejemplo: Urbano Macro, Urbano Micro)
scenario_str = scenario; % Nombre del escenario como cadena para la estructura de resultados

% Estructura para almacenar los datos brutos de SE para el peor usuario (NLOS), por celda
% Esto será L (celdas) x nbrOfSetups (configuraciones) para cada precodificador
results_worst_user_raw = struct();
results_best_users_raw = struct();  
results_all_users_raw = struct();   

%% --- BUCLE PRINCIPAL DE SIMULACIÓN --- %%
for freq_idx = 1:length(frequencies) % Bucle sobre diferentes frecuencias (2,8,15,28)
    frequency = frequencies(freq_idx);
    freq_str = ['f_', strrep(num2str(frequency/1e9), '.', 'p')];
    
    for los_idx = 1:length(LOS_settings) % Bucle sobre diferentes configuraciones LOS (solo LoS=3)
        LoS_setting = LOS_settings(los_idx);
        
        for ris_idx = 1:length(ris_options) % Bucle sobre escenarios con/sin RIS (solo RIS)
            ris_flag = ris_options{ris_idx};
            isRIS = strcmp(ris_flag, 'RIS');
            
            fprintf('Simulando para %s en f = %.0f GHz (LoS=%d %s RIS)\n', ...
                scenario, frequency/1e9, LoS_setting, ris_flag);

            % Inicializar almacenamiento para datos brutos de SE del peor usuario para la combinación freq/LoS/RIS actual
            for prec_idx = 1:length(precoders)
                prec = precoders{prec_idx};
                % Inicializar con ceros: L filas (una por cada celda), nbrOfSetups columnas
                results_worst_user_raw.(freq_str).(ris_flag).(prec) = zeros(L, nbrOfSetups);
                results_best_users_raw.(freq_str).(ris_flag).(prec) = cell(L, nbrOfSetups);
                results_all_users_raw.(freq_str).(ris_flag).(prec) = cell(L, nbrOfSetups);
            end
            
            for n_setup = 1:nbrOfSetups % Bucle sobre diferentes configuraciones de red
                rng(n_setup); % Fijar semilla para la reproducibilidad de cada configuración
                
                if isRIS
                    % --- SIMULACIÓN CON RIS ---
                    [R_BS_UE, HMean_BS_UE, GdB_BS_UE, Kappa_BS_UE, P_LOS_BS_UE, ...
                     R_UE_RIS, HMean_UE_RIS, GdB_UE_RIS, Kappa_UE_RIS, P_LOS_UE_RIS, ...
                     R_BS_RIS_BSAnt, R_BS_RIS_RISel, HMean_BS_RIS, GdB_BS_RIS, Kappa_BS_RIS, P_LOS_BS_RIS, ...
                     worstIdx, bestIdx] = functionExampleSetupRIS(L,K,Mmax,N_ris,ASDdeg,scenario,frequency,LoS_setting,n_setup);
                    
                    GaindB_BS_UE = GdB_BS_UE - noiseVariancedBm;
                    [~,HMean_BS_UE,H_BS_UE_DL_raw,~] = functionChannelGeneration(R_BS_UE, HMean_BS_UE, GaindB_BS_UE, Kappa_BS_UE, P_LOS_BS_UE, K, L, Mmax, nbrOfRealizations);
                    
                    GaindB_BS_RIS = GdB_BS_RIS - noiseVariancedBm;
                    [~,~,~,H_BS_RIS_DL_raw,~] = functionChannelGeneration_BS_RIS(...
                        R_BS_RIS_RISel, ...
                        R_BS_RIS_BSAnt, ...
                        HMean_BS_RIS, ...
                        GaindB_BS_RIS, ...
                        Kappa_BS_RIS, ...
                        P_LOS_BS_RIS, ...
                        N_ris, ...
                        L, ...
                        M, ...
                        nbrOfRealizations);

                    GaindB_UE_RIS = GdB_UE_RIS - noiseVariancedBm;
                    [~,~,H_RIS_UE_DL_raw,~] = functionChannelGeneration(R_UE_RIS, HMean_UE_RIS, GaindB_UE_RIS, Kappa_UE_RIS, P_LOS_UE_RIS, K, L, N_ris, nbrOfRealizations);
                    
                    Theta = calculate_theta(H_BS_UE_DL_raw, H_BS_RIS_DL_raw, H_RIS_UE_DL_raw, N_ris, nbrOfRealizations, L, p, worstIdx);
                    [Heq, Req] = calculate_heq_and_Req(HMean_BS_UE, H_BS_UE_DL_raw, H_BS_RIS_DL_raw, H_RIS_UE_DL_raw, Theta, nbrOfRealizations, Mmax, K, L);
                    
                    [Hhat, C] = functionChannelEstimateMMSE(Req, HMean_BS_UE, Heq, nbrOfRealizations, Mmax, K, L, p, f, tau_p);
                end
                
                [SE_MR, SE_RZF, SE_MMMSE] = functionComputeSE_UL(Hhat, C, Req, tau_c, tau_p, nbrOfRealizations, Mmax, K, L, p, noiseVarianceW);
                
                % Almacenar resultados de SE para el peor usuario por celda directamente en results_worst_user_raw
                for l_cell = 1:L % Bucle a través de cada celda
                    % worstIdx(l_cell) da el índice del usuario peor en la celda l_cell
                    results_worst_user_raw.(freq_str).(ris_flag).MR(l_cell, n_setup) = SE_MR(worstIdx(l_cell), l_cell);
                    results_worst_user_raw.(freq_str).(ris_flag).MMSE(l_cell, n_setup) = SE_MMMSE(worstIdx(l_cell), l_cell);
                    %results_worst_user_raw.(freq_str).(ris_flag).RZF(l_cell, n_setup) = SE_RZF(worstIdx(l_cell), l_cell);
                    
                    best_users = setdiff(1:K, worstIdx(l_cell)); % Los 9 mejores en cada celda
                    results_best_users_raw.(freq_str).(ris_flag).MR{l_cell, n_setup} = SE_MR(best_users, l_cell);
                    results_best_users_raw.(freq_str).(ris_flag).MMSE{l_cell, n_setup} = SE_MMMSE(best_users, l_cell);
                    %results_best_users_raw.(freq_str).(ris_flag).RZF{l_cell, n_setup} = SE_RZF(best_users, l_cell);
                
                    % Almacenar SE para todos los usuarios combinados (LOS & NLOS)
                    results_all_users_raw.(freq_str).(ris_flag).MR{l_cell, n_setup} = SE_MR(:, l_cell);
                    results_all_users_raw.(freq_str).(ris_flag).MMSE{l_cell, n_setup} = SE_MMMSE(:, l_cell);
                    %results_all_users_raw.(freq_str).(ris_flag).RZF{l_cell, n_setup} = SE_RZF(:, l_cell);
                end

                fprintf('Configuración %d (%s, %.0f GHz, LoS=%d %s) Completada.\n', ...
                    n_setup, scenario, frequency/1e9, LoS_setting, ris_flag);
            end
        end
    end
end


%% --- FIGURAS --- %%
colors_prec = containers.Map({'MMMSE','MR'}, {'k','r'});

for freq_idx = 1:length(frequencies)
    plot_frequency = frequencies(freq_idx);
    freq_str_plot = ['f_', strrep(num2str(plot_frequency/1e9), '.', 'p')];
    scenario_plot = scenario;
    ris_flag_plot = ris_options{1};  % Solo "RIS"
    plot_LoS_setting = LOS_settings(1); % Solo LoS=3

    % ---------- Peores Usuarios ----------
    figure('Name', ['CDF: Peores Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']); hold on; box on; grid on;
    for prec = precoders
        data = results_worst_user_raw.(freq_str_plot).(ris_flag_plot).(char(prec));
        data_vec = data(:);
        plot(sort(data_vec), linspace(0,1,length(data_vec))', ...
            [colors_prec(char(prec)) '-'], 'DisplayName', [char(prec) ' (Peor)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - Peores Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']); 
    legend('Location','SouthEast');
    savefig(['CDF_WorstUsers_' freq_str_plot '.fig']);

    % ---------- Mejores Usuarios ----------
    figure('Name', ['CDF: Mejores Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']); hold on; box on; grid on;
    for prec = precoders
        data_all = [];
        for l=1:L
            for s=1:nbrOfSetups
                data_all = [data_all; results_best_users_raw.(freq_str_plot).(ris_flag_plot).(char(prec)){l,s}(:)];
            end    
        end
        plot(sort(data_all), linspace(0,1,length(data_all))', ...
            [colors_prec(char(prec)) '--'], 'DisplayName', [char(prec) ' (Mejor)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - Mejores Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']);
    legend('Location','SouthEast');
    savefig(['CDF_BestUsers_' freq_str_plot '.fig']);

    % ---------- Todos los Usuarios ----------
    figure('Name', ['CDF: Todos los Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']); hold on; box on; grid on;
    for prec = precoders
        data_all = [];
        for l=1:L
            for s=1:nbrOfSetups
                data_all = [data_all; results_all_users_raw.(freq_str_plot).(ris_flag_plot).(char(prec)){l,s}(:)];
            end 
        end
        plot(sort(data_all), linspace(0,1,length(data_all))', ...
            [colors_prec(char(prec)) ':'], 'DisplayName', [char(prec) ' (Todos)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - Todos los Usuarios @ ', num2str(plot_frequency/1e9), ' GHz']);
    legend('Location','SouthEast');
    savefig(['CDF_AllUsers_' freq_str_plot '.fig']);

    % Guardar datos para cada frecuencia
    save(['SE_WorstUser_RIS_data_' freq_str_plot '.mat'], 'results_worst_user_raw');
    save(['SE_BestUsers_RIS_data_' freq_str_plot '.mat'], 'results_best_users_raw');
    save(['SE_AllUsers_RIS_data_' freq_str_plot '.mat'], 'results_all_users_raw');
end

fprintf("Todas las figuras y datos fueron guardados correctamente por frecuencia.\n");
