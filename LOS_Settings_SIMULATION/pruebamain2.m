% Vaciar el espacio de trabajo y cerrar las figuras
close all;
clear;
clc
% Número de BSs
L = 4;
% Número de UEs por BS
K = 10;
% Definir el rango de antenas de la BS
Mrange = 100;
% Extraer el número máximo de antenas de la BS
Mmax = max(Mrange);
% Definir el rango de factores de reutilización de pilotos (no utilizado en este código)
fRange = 1;
% Seleccionar el número de configuraciones con ubicaciones aleatorias de UEs
nbrOfSetups = 100;
% Seleccionar el número de realizaciones de canal por configuración
nbrOfRealizations = 100;
% Ancho de banda de comunicación
B = 20e6;
% Potencia total de transmisión uplink por UE (W)
p = 0.1; % Modificado a 0.1 W para consistencia con SE UL
% Figura de ruido en la BS (en dB)
noiseFigure = 7;
% Calcular la potencia de ruido
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
noiseVarianceW = db2pow(noiseVariancedBm - 30); % Convertir a Watts
% Seleccionar la longitud del bloque de coherencia
tau_c = 200;
% Desviación estándar angular por trayectoria en el modelo de dispersión
% local (en grados)
ASDdeg = 10;
% Parámetro de control de potencia delta (en dB)
deltadB = 10;
% Parámetro de Reutilización de Pilotos
f = 1;
% Longitud del Piloto
tau_p = K;
% Factor prelogaritmo asumiendo solo transmisión UL (Uplink)
prelogFactor = (tau_c - tau_p) / tau_c;
% --- Nuevas variables para selección ---
combinationVector = 'M-MMSE'; % Por defecto (no utilizado directamente)
channelType = 'Ricean'; % Por defecto (implícito en el modelo)
% Definir las frecuencias a simular
frequencies = [2e9, 8e9, 15e9, 28e9]; % en Hz
% Prepara una estructura para guardar los resultados
results = struct();
% Bucle sobre las frecuencias
for freq_idx = 1:length(frequencies)
    frequency = frequencies(freq_idx);
    freq_str = ['f_', strrep(num2str(frequency/1e9), '.', 'p')]; % Añadir el prefijo 'f_'
    % Bucle sobre los escenarios (UMa y UMi)
    for scenario_idx = 1:2
        if scenario_idx == 1
            scenario = 'UMa';
        else
            scenario = 'UMi';
        end
        scenario_str = scenario;
        disp(['Simulando para ' scenario ' a f = ' num2str(frequency/1e9) ' GHz']);
        userSE_MR_NLOS_all = zeros(K*L, nbrOfSetups);
        userSE_M_MMSE_NLOS_all = zeros(K*L, nbrOfSetups);
        userSE_RZF_NLOS_all = zeros(K*L, nbrOfSetups);
        userSE_MR_LOS_all = zeros(K*L, nbrOfSetups);
        userSE_M_MMSE_LOS_all = zeros(K*L, nbrOfSetups);
        userSE_RZF_LOS_all = zeros(K*L, nbrOfSetups);
        userSE_MR_mix_all = zeros(K*L, nbrOfSetups);
        userSE_M_MMSE_mix_all = zeros(K*L, nbrOfSetups);
        userSE_RZF_mix_all = zeros(K*L, nbrOfSetups);
        userSE_MR_case4_all = zeros(K*L, nbrOfSetups);
        userSE_M_MMSE_case4_all = zeros(K*L, nbrOfSetups);
        userSE_RZF_case4_all = zeros(K*L, nbrOfSetups);


        % Bucle sobre las configuraciones
        for n = 1:nbrOfSetups
            % --- Simulación para NLOS (LoS = 0) ---
            LoS = 0;
            [RNormalized_NLOS, HMeanNormalized_NLOS, channelGaindB_NLOS, ricianFactor_NLOS, probLOS_NLOS] = ...
                functionExampleSetup(L, K, Mmax, ASDdeg, scenario, frequency, LoS, n); %n es el índice de la simulación
            channelGainUL_NLOS = channelGaindB_NLOS - noiseVariancedBm;
            userSE_MR_NLOS_setup = zeros(K*L, 1);
            userSE_M_MMSE_NLOS_setup = zeros(K*L, 1);
            userSE_RZF_NLOS_setup = zeros(K*L, 1);
            for m = 1:length(Mrange)
                [R_UL_NLOS, HMean_UL_NLOS, H_UL_NLOS, H_UL_Rayleigh_NLOS] = functionChannelGeneration(RNormalized_NLOS(1:Mrange(m), 1:Mrange(m), :, :, :), ...
                    HMeanNormalized_NLOS(1:Mrange(m), :, :, :), channelGainUL_NLOS, ricianFactor_NLOS, probLOS_NLOS, K, L, Mrange(m), nbrOfRealizations); % Multiplicar ricianFactor por la probabilidad de Rayleigh (1-probLOS)
                [Hhat_MMSE_UL_NLOS, C_MMSE_UL_NLOS] = functionChannelEstimateMMSE(R_UL_NLOS, HMean_UL_NLOS, H_UL_NLOS, nbrOfRealizations, Mrange(m), K, L, p, f, tau_p);
                [SE_MR_matrix_NLOS, SE_RZF_matrix_NLOS, SE_MMMSE_matrix_NLOS] = functionComputeSE_UL(Hhat_MMSE_UL_NLOS, C_MMSE_UL_NLOS, R_UL_NLOS, tau_c, tau_p, nbrOfRealizations, Mrange(m), K, L, p);
                userSE_MR_NLOS_setup = SE_MR_matrix_NLOS(:);
                userSE_M_MMSE_NLOS_setup = SE_MMMSE_matrix_NLOS(:);
                userSE_RZF_NLOS_setup = SE_RZF_matrix_NLOS(:);
                disp(['Configuración ' num2str(n) ' (' scenario ', ' num2str(frequency/1e9) ' GHz, NLOS) para ' num2str(Mmax) ' antenas']);
                clear R_UL_NLOS HMean_UL_NLOS Hhat_MMSE_UL_NLOS C_MMSE_UL_NLOS
            end
            userSE_MR_NLOS_all(:, n) = userSE_MR_NLOS_setup;
            userSE_M_MMSE_NLOS_all(:, n) = userSE_M_MMSE_NLOS_setup;
            userSE_RZF_NLOS_all(:, n) = userSE_RZF_NLOS_setup;
            clear RNormalized_NLOS HMeanNormalized_NLOS channelGaindB_NLOS ricianFactor_NLOS probLOS_NLOS


            % --- Simulación para LOS (LoS = 1) ---
            LoS = 1;
            [RNormalized_LOS, HMeanNormalized_LOS, channelGaindB_LOS, ricianFactor_LOS, probLOS_LOS] = ...
                functionExampleSetup(L, K, Mmax, ASDdeg, scenario, frequency, LoS, n); % La probabilidad de LOS se maneja dentro
            channelGainUL_LOS = channelGaindB_LOS - noiseVariancedBm;
            userSE_MR_LOS_setup = zeros(K*L, 1);
            userSE_M_MMSE_LOS_setup = zeros(K*L, 1);
            userSE_RZF_LOS_setup = zeros(K*L, 1);
            for m = 1:length(Mrange)
                [R_UL_LOS, HMean_UL_LOS, H_UL_LOS, H_UL_Rayleigh_LOS] = functionChannelGeneration(RNormalized_LOS(1:Mrange(m), 1:Mrange(m), :, :, :), ...
                    HMeanNormalized_LOS(1:Mrange(m), :, :, :), channelGainUL_LOS, ricianFactor_LOS, probLOS_LOS, K, L, Mrange(m), nbrOfRealizations); % Multiplicar ricianFactor por la probabilidad de Ricean (probLOS)
                [Hhat_MMSE_UL_LOS, C_MMSE_UL_LOS] = functionChannelEstimateMMSE(R_UL_LOS, HMean_UL_LOS, H_UL_LOS, nbrOfRealizations, Mrange(m), K, L, p, f, tau_p);
                [SE_MR_matrix_LOS, SE_RZF_matrix_LOS, SE_MMMSE_matrix_LOS] = functionComputeSE_UL(Hhat_MMSE_UL_LOS, C_MMSE_UL_LOS, R_UL_LOS, tau_c, tau_p, nbrOfRealizations, Mrange(m), K, L, p);
                userSE_MR_LOS_setup = SE_MR_matrix_LOS(:);
                userSE_M_MMSE_LOS_setup = SE_MMMSE_matrix_LOS(:);
                userSE_RZF_LOS_setup = SE_RZF_matrix_LOS(:);
                disp(['Configuración ' num2str(n) ' (' scenario ', ' num2str(frequency/1e9) ' GHz, LOS) para ' num2str(Mmax) ' antenas']);
                clear R_UL_LOS HMean_UL_LOS Hhat_MMSE_UL_LOS C_MMSE_UL_LOS
            end
            userSE_MR_LOS_all(:, n) = userSE_MR_LOS_setup;
            userSE_M_MMSE_LOS_all(:, n) = userSE_M_MMSE_LOS_setup;
            userSE_RZF_LOS_all(:, n) = userSE_RZF_LOS_setup;
            clear RNormalized_LOS HMeanNormalized_LOS channelGaindB_LOS ricianFactor_LOS probLOS_LOS


            % --- Simulación Mixta (LoS = 2) ---
            LoS = 2;
            [RNormalized_mix, HMeanNormalized_mix, channelGaindB_mix, ricianFactor_mix, probLOS_mix] = ...
                functionExampleSetup(L, K, Mmax, ASDdeg, scenario, frequency, LoS, n); % La probabilidad de LOS se maneja dentro
            channelGainUL_mix = channelGaindB_mix - noiseVariancedBm;
            userSE_MR_mix_setup = zeros(K*L, 1);
            userSE_M_MMSE_mix_setup = zeros(K*L, 1);
            userSE_RZF_mix_setup = zeros(K*L, 1);
            for m = 1:length(Mrange)
                [R_UL_mix, HMean_UL_mix, H_UL_mix, H_UL_Rayleigh_mix] = functionChannelGeneration(RNormalized_mix(1:Mrange(m), 1:Mrange(m), :, :, :), ...
                    HMeanNormalized_mix(1:Mrange(m), :, :, :), channelGainUL_mix, ricianFactor_mix, probLOS_mix, K, L, Mrange(m), nbrOfRealizations); % Multiplicar ricianFactor por la probabilidad de Ricean (probLOS)
                [Hhat_MMSE_UL_mix, C_MMSE_UL_mix] = functionChannelEstimateMMSE(R_UL_mix, HMean_UL_mix, H_UL_mix, nbrOfRealizations, Mrange(m), K, L, p, f, tau_p);
                [SE_MR_matrix_mix, SE_RZF_matrix_mix, SE_MMMSE_matrix_mix] = functionComputeSE_UL(Hhat_MMSE_UL_mix, C_MMSE_UL_mix, R_UL_mix, tau_c, tau_p, nbrOfRealizations, Mrange(m), K, L, p);
                userSE_MR_mix_setup = SE_MR_matrix_mix(:);
                userSE_M_MMSE_mix_setup = SE_MMMSE_matrix_mix(:);
                userSE_RZF_mix_setup = SE_RZF_matrix_mix(:);
                disp(['Configuración ' num2str(n) ' (' scenario ', ' num2str(frequency/1e9) ' GHz, Mixto) para ' num2str(Mmax) ' antenas']);
                clear R_UL_mix HMean_UL_mix Hhat_MMSE_UL_mix C_MMSE_UL_mix
            end
            userSE_MR_mix_all(:, n) = userSE_MR_mix_setup;
            userSE_M_MMSE_mix_all(:, n) = userSE_M_MMSE_mix_setup;
            userSE_RZF_mix_all(:, n) = userSE_RZF_mix_setup;
            clear RNormalized_mix HMeanNormalized_mix channelGaindB_mix ricianFactor_mix probLOS_mix

            
            % --- Simulación para Caso 4 (LoS = 3) ---
            LoS = 3;
            [RNormalized_case4, HMeanNormalized_case4, channelGaindB_case4, ricianFactor_case4, probLOS_case4] = ...
                functionExampleSetup(L, K, Mmax, ASDdeg, scenario, frequency, LoS, n);
            channelGainUL_case4 = channelGaindB_case4 - noiseVariancedBm;
            userSE_MR_case4_setup = zeros(K*L, 1);
            userSE_M_MMSE_case4_setup = zeros(K*L, 1);
            userSE_RZF_case4_setup = zeros(K*L, 1);
            for m = 1:length(Mrange)
                [R_UL_case4, HMean_UL_case4, H_UL_case4, H_UL_Rayleigh_case4] = functionChannelGeneration(RNormalized_case4(1:Mrange(m), 1:Mrange(m), :, :, :), ...
                    HMeanNormalized_case4(1:Mrange(m), :, :, :), channelGainUL_case4, ricianFactor_case4, probLOS_case4, K, L, Mrange(m), nbrOfRealizations);
                [Hhat_MMSE_UL_case4, C_MMSE_UL_case4] = functionChannelEstimateMMSE(R_UL_case4, HMean_UL_case4, H_UL_case4, nbrOfRealizations, Mrange(m), K, L, p, f, tau_p);
                [SE_MR_matrix_case4, SE_RZF_matrix_case4, SE_MMMSE_matrix_case4] = functionComputeSE_UL(Hhat_MMSE_UL_case4, C_MMSE_UL_case4, R_UL_case4, tau_c, tau_p, nbrOfRealizations,Mrange(m), K, L, p);
                userSE_MR_case4_setup = SE_MR_matrix_case4(:);
                userSE_M_MMSE_case4_setup = SE_MMMSE_matrix_case4(:);
                userSE_RZF_case4_setup = SE_RZF_matrix_case4(:);
                disp(['Configuración ' num2str(n) ' (' scenario ', ' num2str(frequency/1e9) ' GHz, Caso 4) para ' num2str(Mmax) ' antenas']);
                clear R_UL_case4 HMean_UL_case4 Hhat_MMSE_UL_case4 C_MMSE_UL_case4
            end
            userSE_MR_case4_all(:, n) = userSE_MR_case4_setup;
            userSE_M_MMSE_case4_all(:, n) = userSE_M_MMSE_case4_setup;
            userSE_RZF_case4_all(:, n) = userSE_RZF_case4_setup;

            % Guardar los resultados para esta configuración
            results.(scenario_str).(freq_str).MR_NLOS = userSE_MR_NLOS_all(:, n);
            results.(scenario_str).(freq_str).MMSE_NLOS = userSE_M_MMSE_NLOS_all(:, n);
            results.(scenario_str).(freq_str).RZF_NLOS = userSE_RZF_NLOS_all(:, n);
            results.(scenario_str).(freq_str).MR_LOS = userSE_MR_LOS_all(:, n);
            results.(scenario_str).(freq_str).MMSE_LOS = userSE_M_MMSE_LOS_all(:, n);
            results.(scenario_str).(freq_str).RZF_LOS = userSE_RZF_LOS_all(:, n);
            results.(scenario_str).(freq_str).MR_mix = userSE_MR_mix_all(:, n);
            results.(scenario_str).(freq_str).MMSE_mix = userSE_M_MMSE_mix_all(:, n);
            results.(scenario_str).(freq_str).RZF_mix = userSE_RZF_mix_all(:, n);
            results.(scenario_str).(freq_str).MR_case4 = userSE_MR_case4_all(:, n);
            results.(scenario_str).(freq_str).MMSE_case4 = userSE_M_MMSE_case4_all(:, n);
            results.(scenario_str).(freq_str).RZF_case4 = userSE_RZF_case4_all(:, n);

            clear RNormalized_NLOS HMeanNormalized_NLOS channelGaindB_NLOS ricianFactor_NLOS probLOS_NLOS
            clear RNormalized_LOS HMeanNormalized_LOS channelGaindB_LOS ricianFactor_LOS probLOS_LOS
            clear RNormalized_mix HMeanNormalized_mix channelGaindB_mix ricianFactor_mix probLOS_mix
            clear RNormalized_case4 HMeanNormalized_case4 channelGaindB_case4 ricianFactor_case4 probLOS_case4
        end
        % Promediar los resultados sobre las configuraciones
        results.(scenario_str).(freq_str).MR_NLOS = mean(userSE_MR_NLOS_all, 2);
        results.(scenario_str).(freq_str).MMSE_NLOS = mean(userSE_M_MMSE_NLOS_all, 2);
        results.(scenario_str).(freq_str).RZF_NLOS = mean(userSE_RZF_NLOS_all, 2);
        results.(scenario_str).(freq_str).MR_LOS = mean(userSE_MR_LOS_all, 2);
        results.(scenario_str).(freq_str).MMSE_LOS = mean(userSE_M_MMSE_LOS_all, 2);
        results.(scenario_str).(freq_str).RZF_LOS = mean(userSE_RZF_LOS_all, 2);
        results.(scenario_str).(freq_str).MR_mix = mean(userSE_MR_mix_all, 2);
        results.(scenario_str).(freq_str).MMSE_mix = mean(userSE_M_MMSE_mix_all, 2);
        results.(scenario_str).(freq_str).RZF_mix = mean(userSE_RZF_mix_all, 2);
        results.(scenario_str).(freq_str).MR_case4 = mean(userSE_MR_case4_all, 2);
        results.(scenario_str).(freq_str).MMSE_case4 = mean(userSE_M_MMSE_case4_all, 2);
        results.(scenario_str).(freq_str).RZF_case4 = mean(userSE_RZF_case4_all, 2);
    end
end
%% Graficar las figuras
CDFnumbers = linspace(0, 1, K * L); % Ahora la CDF es sobre los usuarios promedio
figureOffset = 0;
frequencies_plot = [2e9, 8e9, 15e9, 28e9];
scenarios_plot = {'UMa', 'UMi'};
colors = {'r', 'k', 'b', 'g'}; % Añadido verde para el caso 4
linestyles = {'--', '-', ':', '-.'}; % -- para NLOS, - para LOS, : para Mixto, -. para Caso 4
receiver_schemes = {'MR', 'MMSE', 'RZF'};
for i = 1:length(scenarios_plot)
    scenario = scenarios_plot(i);
    for j = 1:length(frequencies_plot)
        frequency = frequencies_plot(j);
        freq_str = ['f_', strrep(num2str(frequency/1e9), '.', 'p')];
        figure(i * length(frequencies_plot) - (length(frequencies_plot) - j) + figureOffset);
        hold on;
        box on;
        for k = 1:length(receiver_schemes)
            scheme = receiver_schemes{k};
            % Traza NLOS, LOS, Mixto y Caso 4 en la misma figura con diferentes estilos de línea
            plot(sort(results.(scenario{1}).(freq_str).([scheme '_NLOS'])), CDFnumbers, 'Color', colors{k}, 'LineStyle', linestyles{1}, 'DisplayName', [scheme ' (NLOS)']);
            plot(sort(results.(scenario{1}).(freq_str).([scheme '_LOS'])), CDFnumbers, 'Color', colors{k}, 'LineStyle', linestyles{2}, 'DisplayName', [scheme ' (LOS)']);
            plot(sort(results.(scenario{1}).(freq_str).([scheme '_mix'])), CDFnumbers, 'Color', colors{k}, 'LineStyle', linestyles{3}, 'DisplayName', [scheme ' (Mixto)']);
            plot(sort(results.(scenario{1}).(freq_str).([scheme 'LOS 9x1 NLOS'])), CDFnumbers, 'Color', colors{k}, 'LineStyle', linestyles{4}, 'DisplayName', [scheme ' (LoS 9x1 NLOS)']);
        end
        xlabel('SE por UE [bit/s/Hz]');
        ylabel('Función de Distribución Acumulada (CDF)');
        %         title([scenario{1} ' a f = ' num2str(frequency/1e9) ' GHz']); % Eliminar el título
        legend('Location', 'SouthEast');
        grid on;
        % Guardar la figura
        filename = sprintf('%s_f_%dGHz_CDF.png', scenario{1}, round(frequency/1e9));
        saveas(gcf, filename); % Guarda la figura como un archivo .png
    end
end
