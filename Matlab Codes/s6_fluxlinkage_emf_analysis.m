%% fluxlinkage_emf_analysis.m — λ & EMF KPIs + spectra (with corrected E1/THD)
% Post-process: flux linkages and EMF from FEMM_PostResults56.mat
% Computes:
%   - Phase & line-line RMS EMF
%   - Fundamental E1 and THD (phase A) via projection
%   - EMF spectrum
%   - Saves KPIs & spectra and generates plots

% Run in the folder containing FEMM_PostResults56.mat

% clear; 
clc;
% close all;

%% ---------------------- User Parameters ----------------------
inMat  = 'FEMM_PostResults.mat';   % Input from s4_master_sweep_and_postproc.m
outMat = 'emf_analysis.mat';       % Output .mat file
outCsv = 'emf_summary.csv';        % Output .csv summary
Npk    = 5;                          % Number of largest spectral peaks to label

%% ---------------------- Load Results -------------------------
S = load(inMat);
R = S.results;

ang_m = R.angle_mech(:);             % mechanical angle (deg)

% Electrical angle: use saved value if present, otherwise reconstruct p*theta_mech
if isfield(R,'angle_elec')
    ang_e = R.angle_elec(:);
else
    ang_e = R.meta.p * ang_m(:);
end

p           = R.meta.p;
omega_mech  = R.meta.omega_mech;
stepMechDeg = R.meta.stepMechDeg;
dtheta      = deg2rad(stepMechDeg);  % mechanical step [rad]

lamA = R.lambda.A(:);
lamB = R.lambda.B(:);
lamC = R.lambda.C(:);

%% ---------------------- EMF Compute (if missing) ------------

if isfield(R,'emf') && isfield(R.emf,'A')
    eA = R.emf.A(:);
    eB = R.emf.B(:);
    eC = R.emf.C(:);
else
    dlamA = (circshift(lamA,-1) - circshift(lamA,1)) / (2*dtheta);
    dlamB = (circshift(lamB,-1) - circshift(lamB,1)) / (2*dtheta);
    dlamC = (circshift(lamC,-1) - circshift(lamC,1)) / (2*dtheta);
    eA = -dlamA * omega_mech;
    eB = -dlamB * omega_mech;
    eC = -dlamC * omega_mech;
end

%% ---------------------- KPIs (phase & line-line RMS) --------
Erms_A = rms(eA);
Erms_B = rms(eB);
Erms_C = rms(eC);

eAB = eA - eB;
eBC = eB - eC;
eCA = eC - eA;

Erms_AB = rms(eAB);
Erms_BC = rms(eBC);
Erms_CA = rms(eCA);

%% -------- Fundamental (E1) & THD via projection -------------
theta_e_rad = deg2rad( mod(ang_e - ang_e(1), 360) );   % 0..2π over one elec period
eA0 = eA - mean(eA);                                   % remove DC

s1 = sin(theta_e_rad);
c1 = cos(theta_e_rad);

a = (2/numel(eA0)) * sum(eA0 .* s1);
b = (2/numel(eA0)) * sum(eA0 .* c1);

E1_amp = hypot(a,b);                % fundamental peak amplitude
E1_rms = E1_amp / sqrt(2);          % fundamental RMS
e1     = a*s1 + b*c1;               % reconstructed fundamental (peak)
harm   = eA0 - e1;                  % all harmonics (AC) except order-1
THD_A  = 100 * ( rms(harm) / E1_rms );

%% ---------------------- Print Summary -----------------------
fprintf('\nBack-EMF KPIs\n');
fprintf('Phase RMS (V):      A=%.3f  B=%.3f  C=%.3f\n', Erms_A,Erms_B,Erms_C);
fprintf('Line-line RMS (V): AB=%.3f  BC=%.3f  CA=%.3f\n', Erms_AB,Erms_BC,Erms_CA);
fprintf('Corrected fundamental (phase A): E1_RMS=%.3f V, THD=%.2f %%\n', E1_rms, THD_A);

%% ---------------------- Spectrum (Phase A) ------------------
x = eA0;                            % AC component (no DC)
N = numel(x);
w = hann(N);
X = fft(x .* w);
W = sum(w) / 2;                     % amplitude norm (for relative use)
magA   = abs(X(1:floor(N/2)+1)) / (W * N);   % single-sided amplitude-ish
orders = (0:floor(N/2)).';                   % electrical order index (0..)

%% ---------------------- Save Results ------------------------
emf_kpis = struct( ...
    'Erms_phase', [Erms_A,Erms_B,Erms_C], ...
    'Erms_ll',    [Erms_AB,Erms_BC,Erms_CA], ...
    'E1_rms_A',   E1_rms, ...
    'THD_A_pct',  THD_A);

emf_spectrum = struct('orders',orders,'magA',magA);

save(outMat, ...
     'ang_m','ang_e', ...
     'lamA','lamB','lamC', ...
     'eA','eB','eC','eAB','eBC','eCA', ...
     'emf_kpis','emf_spectrum');

fid = fopen(outCsv,'w');
fprintf(fid,'ErmsA,ErmsB,ErmsC,ErmsAB,ErmsBC,ErmsCA,E1rmsA,THD_A_pct\n');
fprintf(fid,'%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f\n', ...
        Erms_A,Erms_B,Erms_C,Erms_AB,Erms_BC,Erms_CA,E1_rms,THD_A);
fclose(fid);

%% ---------------------- Plots -------------------------------
figure;
plot(ang_m, [lamA lamB lamC], 'LineWidth', 1); grid on;
xlabel('\theta_{mech} (deg)'); ylabel('\lambda (Wb-turn)');
title('Flux Linkage per Phase'); legend('A','B','C');

figure;
plot(ang_m, [eA eB eC], 'LineWidth', 1); grid on;
xlabel('\theta_{mech} (deg)'); ylabel('e (V)');
title('Back-EMF per Phase'); legend('e_A','e_B','e_C');

figure;
plot(ang_m, [eAB eBC eCA], 'LineWidth', 1); grid on;
xlabel('\theta_{mech} (deg)'); ylabel('e_{LL} (V)');
title('Line-to-Line Back-EMF'); legend('e_{AB}','e_{BC}','e_{CA}');

figure;
stem(orders, magA, 'filled'); grid on;
xlabel('Electrical order'); ylabel('Phase-A amplitude (rel.)');
title('Back-EMF Spectrum (Phase A)');

% Label top-Npk peaks
[~, idxTop] = maxk(magA, Npk);
text(orders(idxTop) + 0.8, magA(idxTop), string(orders(idxTop)), ...
     'FontWeight','bold');
