%% main2.m - Script principal para Massive MIMO sin RIS (UMi)
close all; clear; clc;

%% Parámetros de la simulación
L = 4; 
K = 10; 
M = 40; 
Mmax = M;
nbrOfSetups = 100; 
nbrOfRealizations = 100;
B = 20e6; 
p = 0.1; 
noiseFigure = 7;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
noiseVarianceW = db2pow(noiseVariancedBm - 30);
tau_c = 200; 
tau_p = K;
ASDdeg = 10; 
f = 1;
frequencies = [2e9, 8e9, 28e9];
precoders_to_run = {'MR', 'MMMSE', 'RZF'};

scenario = 'UMi';

results_LOS = struct();
results_NLOS = struct();
results_ALL = struct();

colors = containers.Map({'MR', 'MMMSE', 'RZF'}, {'r', 'k', 'b'});

%% Bucle por frecuencia
for freq_idx = 1:length(frequencies)
    frequency = frequencies(freq_idx);
    freq_str = ['f_' strrep(num2str(frequency/1e9), '.', 'p')];
    disp(['Simulando ' scenario ' a ' num2str(frequency/1e9) ' GHz sin RIS']);

    for prec = precoders_to_run
        precoder = prec{1};
        results_LOS.(precoder) = [];
        results_NLOS.(precoder) = [];
        results_ALL.(precoder) = [];
    end

    for n_setup = 1:nbrOfSetups
        LoS_setting = 3;

        [R_BS_UE, HMean_BS_UE, channelGaindB_BS_UE, ricianFactor_BS_UE, probLOS_BS_UE, ...
         worstUserIndex_per_cell, bestUsersIndices_per_cell] = ...
            functionExampleSetup(L, K, Mmax, ASDdeg, scenario, frequency, LoS_setting, n_setup);

        channelGainUL_BS_UE = channelGaindB_BS_UE - noiseVariancedBm;

        [R_UL_BS_UE, HMean_UL_BS_UE, H_UL_BS_UE] = ...
            functionChannelGeneration(R_BS_UE, HMean_BS_UE, channelGainUL_BS_UE, ...
                                      ricianFactor_BS_UE, probLOS_BS_UE, ...
                                      K, L, Mmax, nbrOfRealizations);

        [Hhat, C] = functionChannelEstimateMMSE(R_UL_BS_UE, HMean_UL_BS_UE, H_UL_BS_UE, ...
                                                nbrOfRealizations, Mmax, K, L, p, f, tau_p);

        [SE_MR, SE_RZF, SE_MMMSE] = functionComputeSE_UL(Hhat, C, R_UL_BS_UE, ...
                            tau_c, tau_p, nbrOfRealizations, Mmax, K, L, p, noiseVarianceW);

        prec_dict = struct('MR', SE_MR, 'RZF', SE_RZF, 'MMMSE', SE_MMMSE);

        for prec = precoders_to_run
            precoder = prec{1};
            SE_matrix = prec_dict.(precoder);

            se_nlos = zeros(1, L);
            se_los = [];
            se_all = [];

            for l_cell = 1:L
                se_nlos(l_cell) = SE_matrix(worstUserIndex_per_cell(l_cell), l_cell);
                se_los = [se_los; SE_matrix(bestUsersIndices_per_cell{l_cell}, l_cell)];
                se_all = [se_all; SE_matrix(:, l_cell)];
            end

            results_NLOS.(precoder) = [results_NLOS.(precoder); se_nlos(:)];
            results_LOS.(precoder) = [results_LOS.(precoder); se_los(:)];
            results_ALL.(precoder) = [results_ALL.(precoder); se_all(:)];
        end

        disp(['Setup ' num2str(n_setup) ' completado (f = ' num2str(frequency/1e9) ' GHz)']);
    end

    %% GUARDAR DATOS
    save('SE_NLOS_noRIS_data.mat', 'results_NLOS');
    save('SE_LOS_noRIS_data.mat', 'results_LOS');
    save('SE_ALL_noRIS_data.mat', 'results_ALL');

    %% GRAFICAR CDFs
    figure('Name', 'CDF NLOS'); hold on; box on; grid on;
    for prec = precoders_to_run
        precoder = prec{1};
        data = sort(results_NLOS.(precoder));
        plot(data, linspace(0,1,length(data))', [colors(precoder) '--'], 'DisplayName', [precoder ' (NLOS)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - NLOS Users (Worst) @ ' num2str(frequency/1e9) ' GHz']);
    legend('Location','SouthEast');
    saveas(gcf, ['CDF_NLOS_' freq_str '.png']);
    savefig(['CDF_NLOS_' freq_str '.fig']);

    figure('Name', 'CDF LOS'); hold on; box on; grid on;
    for prec = precoders_to_run
        precoder = prec{1};
        data = sort(results_LOS.(precoder));
        plot(data, linspace(0,1,length(data))', [colors(precoder) '-.'], 'DisplayName', [precoder ' (LOS)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - LOS Users (Best) @ ' num2str(frequency/1e9) ' GHz']);
    legend('Location','SouthEast');
    saveas(gcf, ['CDF_LOS_' freq_str '.png']);
    savefig(['CDF_LOS_' freq_str '.fig']);

    figure('Name', 'CDF ALL'); hold on; box on; grid on;
    for prec = precoders_to_run
        precoder = prec{1};
        data = sort(results_ALL.(precoder));
        plot(data, linspace(0,1,length(data))', [colors(precoder) ':'], 'DisplayName', [precoder ' (All)']);
    end
    xlabel('SE [bit/s/Hz]'); ylabel('CDF');
    title(['CDF - ALL Users @ ' num2str(frequency/1e9) ' GHz']);
    legend('Location','SouthEast');
    saveas(gcf, ['CDF_ALL_' freq_str '.png']);
    savefig(['CDF_ALL_' freq_str '.fig']);
end