%%main_RIS_No_RIS
close all;
clear;
clc;
%% --- USER CONFIGURABLE INPUTS --- %%
frequencies = 2e9;    % Frequencies (Hz) for simulation (fixed to 15 GHz as requested)
LOS_settings = 3;      % Line-of-Sight (LOS) configurations (fixed to 3 as requested)
precoders = ["MR", "MMSE", "RZF"];      % Precoding schemes to simulate
ris_options = "RIS";        % Scenarios to simulate: with RIS or without (fixed to RIS as requested)

%% --- FIXED SIMULATION PARAMETERS --- %%
L = 4;                 % Number of cells (Base Stations - BS)
K = 10;                % Number of users (User Equipments - UE) per cell
M = 40;                % Number of antennas at each BS (M = Mmax)
Mmax = M;              % Alias for maximum number of antennas at the BS
N_ris = 100;            % Number of elements in the RIS
nbrOfSetups = 100;       % Number of network setups to average over
nbrOfRealizations = 100; % Number of channel realizations per setup (for averaging within a setup)
B = 20e6;              % System bandwidth (Hz)
p = 0.1;               % Transmit power per UE (W)
noiseFigure = 7;       % Receiver noise figure (dB)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure; % Noise variance in dBm
noiseVarianceW = db2pow(noiseVariancedBm - 30); % Noise variance in Watts
tau_c = 200;           % Channel coherence length (symbols)
tau_p = K;             % Pilot sequence length (symbols)
f = 1;                 % Factor for noise power calculation (typically 1 for AWGN noise)
ASDdeg = 10;           % Angular Spread of Departures/Arrivals (degrees)
scenario = 'UMi';      % Propagation scenario (e.g. Urban Macro, Urban Micro)
scenario_str = scenario; % Scenario name as string for results structure

% Structure to store the raw SE data for the worst (NLOS) users, per cell
% This will be L (cells) x nbrOfSetups (setups) for each precoder
results_worst_user_raw = struct();
results_best_users_raw = struct();  
results_all_users_raw = struct();   

%% --- MAIN SIMULATION LOOP --- %%
for freq_idx = 1:length(frequencies) % Loop over different frequencies (only 15 GHz)
    frequency = frequencies(freq_idx);
    freq_str = ['f_', strrep(num2str(frequency/1e9), '.', 'p')];
    
    for los_idx = 1:length(LOS_settings) % Loop over different LOS configurations (only LoS=3)
        LoS_setting = LOS_settings(los_idx);
        
        for ris_idx = 1:length(ris_options) % Loop over scenarios with/without RIS (only RIS)
            ris_flag = ris_options{ris_idx};
            isRIS = strcmp(ris_flag, 'RIS');
            
            fprintf('Simulating for %s at f = %.0f GHz (LoS=%d %s RIS)\n', ...
                scenario, frequency/1e9, LoS_setting, ris_flag);

            % Initialize storage for raw worst-user SE data for the current freq/LoS/RIS combo
            for prec_idx = 1:length(precoders)
                prec = precoders{prec_idx};
                % Initialize with zeros: L rows (one for each cell), nbrOfSetups columns
                results_worst_user_raw.(freq_str).(ris_flag).(prec) = zeros(L, nbrOfSetups);
                results_best_users_raw.(freq_str).(ris_flag).(prec) = cell(L, nbrOfSetups);
                results_all_users_raw.(freq_str).(ris_flag).(prec) = cell(L, nbrOfSetups);
            end
            
            for n_setup = 1:nbrOfSetups % Loop over different network setups
                rng(n_setup); % Set seed for reproducibility of each setup
                
                if isRIS
                    % --- SIMULATION WITH RIS ---
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
                
                [SE_MR, SE_RZF, SE_MMSE] = functionComputeSE_UL(Hhat, C, Req, tau_c, tau_p, nbrOfRealizations, Mmax, K, L, p, noiseVarianceW);
                
                % Store SE results for the worst user per cell directly into results_worst_user_raw
                for l_cell = 1:L % Loop through each cell
                    % worstIdx(l_cell) gives the user index for the worst user in cell l_cell
                    results_worst_user_raw.(freq_str).(ris_flag).MR(l_cell, n_setup) = SE_MR(worstIdx(l_cell), l_cell);
                    results_worst_user_raw.(freq_str).(ris_flag).MMSE(l_cell, n_setup) = SE_MMSE(worstIdx(l_cell), l_cell);
                    results_worst_user_raw.(freq_str).(ris_flag).RZF(l_cell, n_setup) = SE_RZF(worstIdx(l_cell), l_cell);
                    
                    best_users = setdiff(1:K, worstIdx(l_cell)); % Los 9 mejores en cada celda
                    results_best_users_raw.(freq_str).(ris_flag).MR{l_cell, n_setup} = SE_MR(best_users, l_cell);
                    results_best_users_raw.(freq_str).(ris_flag).MMSE{l_cell, n_setup} = SE_MMSE(best_users, l_cell);
                    results_best_users_raw.(freq_str).(ris_flag).RZF{l_cell, n_setup} = SE_RZF(best_users, l_cell);
                
                    % Almacenar todos los usuarios (para posible análisis global)
                    results_all_users_raw.(freq_str).(ris_flag).MR{l_cell, n_setup} = SE_MR(:, l_cell);
                    results_all_users_raw.(freq_str).(ris_flag).MMSE{l_cell, n_setup} = SE_MMSE(:, l_cell);
                    results_all_users_raw.(freq_str).(ris_flag).RZF{l_cell, n_setup} = SE_RZF(:, l_cell);
                end

                fprintf('Configuración %d (%s, %.0f GHz, LoS=%d %s) Completada.\n', ...
                    n_setup, scenario, frequency/1e9, LoS_setting, ris_flag);
            end
        end
    end
end

%%Graficar las figuras:
%% --- GRÁFICAS --- %%
colors_prec = containers.Map({'MMSE','RZF','MR'}, {'k','b','r'});

plot_frequency = frequencies(1);
freq_str_plot = ['f_', strrep(num2str(plot_frequency/1e9), '.', 'p')];
scenario_plot = scenario;
ris_flag_plot = ris_options{1};
plot_LoS_setting = LOS_settings(1);

% ---------- Worst Users ----------
figure('Name', 'CDF: Worst Users'); hold on; box on; grid on;
for prec = precoders
    data = results_worst_user_raw.(freq_str_plot).(ris_flag_plot).(char(prec));
    data_vec = data(:);
    plot(sort(data_vec), linspace(0,1,length(data_vec))', ...
        [colors_prec(char(prec)) '-'], 'DisplayName', [char(prec) ' (Worst)']);
end
xlabel('SE [bit/s/Hz]'); ylabel('CDF');
title('CDF - Worst Users'); legend('Location','SouthEast');
savefig('CDF_WorstUsers.fig');

% ---------- Best Users ----------
figure('Name', 'CDF: Best Users'); hold on; box on; grid on;
for prec = precoders
    data_all = [];
    for l=1:L
        for s=1:nbrOfSetups
            data_all = [data_all; results_best_users_raw.(freq_str_plot).(ris_flag_plot).(char(prec)){l,s}(:)];
        end    
    end
    plot(sort(data_all), linspace(0,1,length(data_all))', ...
        [colors_prec(char(prec)) '--'], 'DisplayName', [char(prec) ' (Best)']);
end
xlabel('SE [bit/s/Hz]'); ylabel('CDF');
title('CDF - Best Users'); legend('Location','SouthEast');
savefig('CDF_BestUsers.fig');

% ---------- All Users ----------
figure('Name', 'CDF: All Users'); hold on; box on; grid on;
for prec = precoders
    data_all = [];
    for l=1:L
        for s=1:nbrOfSetups
            data_all = [data_all; results_all_users_raw.(freq_str_plot).(ris_flag_plot).(char(prec)){l,s}(:)];
        end 
    end
    plot(sort(data_all), linspace(0,1,length(data_all))', [colors_prec(char(prec)) ':'], 'DisplayName', [char(prec) ' (All)']);
end
xlabel('SE [bit/s/Hz]'); ylabel('CDF');
title('CDF - All Users'); legend('Location','SouthEast');
savefig('CDF_AllUsers.fig');

%% --- GUARDAR DATOS --- %%
save('SE_WorstUser_RIS_data.mat', 'results_worst_user_raw');
save('SE_BestUsers_RIS_data.mat', 'results_best_users_raw');
save('SE_AllUsers_RIS_data.mat', 'results_all_users_raw');
fprintf("Datos guardados correctamente.\n");