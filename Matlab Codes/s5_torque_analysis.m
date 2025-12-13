%% torque_analysis.m — KPIs + spectrum (with harmonic markers)
% Post-processes FEMM sweep results:
%   - Computes torque ripple KPIs
%   - Computes torque spectrum over one electrical period
%   - Saves .mat and .csv summaries
%   - Plots torque vs mechanical angle and its spectrum

clc;
% clear;   % keep commented so workspace (e.g. design data) isn't destroyed

%% ---------------------- User Parameters ----------------------
inMat  = 'FEMM_PostResults.mat';            % Input results from master_sweep_and_postproc
outMat = 'torque_analysis.mat';             % Output .mat file (KPIs + spectrum)
outCsv = 'torque_analysis_summary.csv';     % Output .csv summary

Npk  = 5;                                     % Number of largest spectral peaks to label

%% ---------------------- Load Results -------------------------
S = load(inMat);
R = S.results;

p      = R.meta.p;                            % pole pairs from master
T      = R.torque(:);                         % torque waveform [N·m]
ang_m  = R.angle_mech(:);                     % mechanical angle [deg]

% Electrical angle: use saved value if present, otherwise reconstruct p*theta_mech
if isfield(R,'angle_elec')
    ang_e = R.angle_elec(:);                  % electrical angle [deg]
else
    ang_e = p * ang_m(:);                     % fallback (one electrical period assumption)
end

%% ---------------------- Torque KPIs -------------------------
Tavg        = mean(T);
Tpp         = max(T) - min(T);
Tripple_pct = 100 * Tpp / max(1e-9, abs(Tavg));          % peak-to-peak ripple vs |Tavg|
Trms_ac     = rms(T - Tavg);                             % AC component RMS
THD_vs_DC_pct = 100 * Trms_ac / max(1e-9, abs(Tavg));    % RMS ripple vs |Tavg|

fprintf('\nTorque KPIs\n');
fprintf('T_avg       = %.6f N·m\n', Tavg);
fprintf('T_p2p       = %.6f N·m\n', Tpp);
fprintf('Ripple p2p  = %.3f %% of |T_avg|\n', Tripple_pct);
fprintf('Ripple RMS  = %.6f N·m (%.3f %% of |T_avg|)\n', Trms_ac, THD_vs_DC_pct);

%% ---------------------- Spectrum ----------------------------
x = T - Tavg;                         % remove DC
N = numel(x);
w = hann(N);                          % Hann window
X = fft(x .* w);
W = sum(w) / 2;
mag    = abs(X(1:floor(N/2)+1)) / (W * N);   % single-sided amplitude-like spectrum
orders = (0:floor(N/2)).';                   % electrical orders (one electrical period assumed)

%% ---------------------- Save Numbers ------------------------
torque_kpis = struct( ...
    'Tavg',                 Tavg, ...
    'Tpp',                  Tpp, ...
    'RipplePct_p2p',        Tripple_pct, ...
    'RippleRMS',            Trms_ac, ...
    'RippleRMS_over_DC_pct',THD_vs_DC_pct);

torque_spectrum = struct('orders',orders,'mag',mag);

save(outMat,'torque_kpis','torque_spectrum','ang_m','ang_e','T');

fid = fopen(outCsv,'w');
fprintf(fid,'Tavg_Nm,Tp2p_Nm,RipplePct_p2p,RippleRMS_Nm,RippleRMS_over_DC_pct\n');
fprintf(fid,'%.9f,%.9f,%.6f,%.9f,%.6f\n',Tavg,Tpp,Tripple_pct,Trms_ac,THD_vs_DC_pct);
fclose(fid);

%% ---------------------- Plots ------------------------------
% Torque vs mechanical angle
figure;
plot(ang_m, T, 'LineWidth', 1); grid on;
xlabel('\theta_{mech} (deg)'); ylabel('Torque (N·m)');
title('Torque vs Mechanical Angle');

% Spectrum with labeled peaks
figure;
stem(orders, mag, 'filled'); grid on;
xlabel('Electrical order'); ylabel('Amplitude (N·m)');
title('Torque Ripple Spectrum'); hold on;

% (optional) your vertical guides using Ns, p etc. can go here
% e.g. slotting orders, 2p, 6p, Ns±2p, etc.

% Label top-N peaks by order (requires Signal Processing Toolbox: findpeaks)
[pk, idx] = findpeaks(mag, 'SortStr','descend', 'NPeaks', Npk);

for i = 1:numel(idx)
    text(orders(idx(i)) + 0.8, pk(i), sprintf('%d', orders(idx(i))), ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', [0.2 0.2 0.2]);
end

hold off;