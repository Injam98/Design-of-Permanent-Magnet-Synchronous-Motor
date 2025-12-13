%% master_sweep_and_postproc.m — torque ripple style stepping + post-proc
% Requires:
%   - s1_geometry_calc.m (analytical design) already run
%   - s2_FEMM_Geometry.m (FEMM model build) already run
%   - s3_torque_vs_thetae.m (to determine theta_e_star_deg)
%
% So that the following exist in the workspace:
%   p, rpm_mech, I_d, I_q, femm_path, GROUP_ROTOR, femm_filename, theta_e_star_deg

%clc; clear; close all;

%% ---------------------- User Inputs ----------------------
% p, rpm_mech come from Script 1

stepElecDeg    = 1;                          % electrical step [deg]

% Use the FEM file generated in Script 2 as base
baseFemFile    = femm_filename;              % original model (from s2_FEMM_Geometry.m)
tmpFemFile     = 'PMSM_TMP.fem';             % working copy (auto-deleted)
outMat         = 'FEMM_PostResults.mat';   % results file

% ---- Current control ----
mode           = 'load';                     % 'load' (torque-producing) or 'open' (cogging/EMF)

Id             = I_d;                        % d-axis current [A] (from Script 1)
Iq             = I_q;                        % q-axis current [A] (from Script 1)

% θe* should be computed by Script 3 (torque_vs_thetae.m)
if ~exist('theta_e_star_deg','var')
    error('theta_e_star_deg not found in workspace. Run s3_torque_vs_thetae.m first.');
end

% ---- Sweep length (mech degrees). For one electrical period: 360/p
rotationMechDeg = 360/p;                     % e.g., 360/p (one elec period). You can set 180, 360, etc.

%% ---------------------- Derived Settings -----------------
stepMechDeg = stepElecDeg / p;               % mechanical step [deg]
nSteps      = round(rotationMechDeg / stepMechDeg);

if abs(rotationMechDeg/stepMechDeg - nSteps) > 1e-12
    error('rotationMechDeg must be a multiple of stepMechDeg.');
end

theta_mech  = (0:nSteps-1) * stepMechDeg;    % mech angle vector [deg], no duplicate end
theta_elec  = theta_mech * p;                % electrical angle [deg]
omega_mech  = rpm_mech * 2*pi/60;            % [rad/s]
dtheta_rad  = deg2rad(stepMechDeg);          % mech step [rad]

fprintf('p=%d | stepElec=%g° | stepMech=%g° | sweep=%g° mech | steps=%d | rpm=%g | θe*=%g°\n', ...
    p, stepElecDeg, stepMechDeg, rotationMechDeg, nSteps, rpm_mech, theta_e_star_deg);

%% ---------------------- FEMM Setup (hidden) ---------------
% femm_path is defined in Script 2 (s2_FEMM_Geometry.m)
addpath(femm_path);
openfemm(1);                                 % run FEMM hidden (background)
opendocument(baseFemFile);                   % open the same model created by Script 2
mi_saveas(tmpFemFile);                       % working copy; original stays untouched
mi_smartmesh(0);                             % quick mesh
fprintf('Working file: %s\n', tmpFemFile);

%% ---------------------- Preallocate -----------------------
T_mech    = zeros(nSteps,1);                      % torque [Nm]
lamA = zeros(nSteps,1);
lamB = lamA;
lamC = lamA;   % flux linkages [Wb-turn]

%% ---------------------- Sweep (single FEMM session) -------
mi_seteditmode('group');
fprintf('\nStarting rotor sweep...\n');

for k = 1:nSteps
    th_mech = theta_mech(k);                 % current mech angle [deg]
    th_elec = theta_e_star_deg + theta_elec(k);  % synchronous electrical angle [deg]

    % --- Set phase currents (either open-circuit or dq->abc following rotor) ---
    if strcmpi(mode,'open')                 % Set 'load' for torque or 'open' for cogging torque
        Ia = 0;
        Ib = 0;
        Ic = 0;
    else
        the = deg2rad(th_elec);             % electrical angle [rad]
        Ia  = sqrt(2) * ( Id*cos(the)               - Iq*sin(the) );
        Ib  = sqrt(2) * ( Id*cos(the - 2*pi/3)      - Iq*sin(the - 2*pi/3) );
        Ic  = sqrt(2) * ( Id*cos(the + 2*pi/3)      - Iq*sin(the + 2*pi/3) );
    end
    % Your convention: ('A',1,current)
    mi_modifycircprop('A',  1, Ia);
    mi_modifycircprop('B',  1, Ib);
    mi_modifycircprop('C',  1, Ic);

    % --- Solve & load (silent) ---
    mi_analyze(false);
    mi_loadsolution;

    % --- Torque from rotor groups ---
    % GROUP_ROTOR is defined in Script 2
    mo_groupselectblock(GROUP_ROTOR);
    for gi = 1:(2*p), mo_groupselectblock(GROUP_ROTOR+gi); end
    T_mech(k) = mo_blockintegral(22);
    mo_clearblock;

    % --- Phase flux linkages λ (3rd value from circuit props) ---
    a = mo_getcircuitproperties('A'); lamA(k) = a(3);
    b = mo_getcircuitproperties('B'); lamB(k) = b(3);
    c = mo_getcircuitproperties('C'); lamC(k) = c(3);

    % --- Status print: step, mech°, elec°, currents, torque ---
    fprintf('Step %3d / %3d | θ_mech=%7.2f° | θ_elec=%7.2f° | Ia=%7.1f Ib=%7.1f Ic=%7.1f | T=%9.3f Nm\n', ...
        k, nSteps, th_mech, th_elec, Ia, Ib, Ic, T_mech(k));

    % --- Advance rotor to next step (skip after last step) ---
    if k < nSteps
        mi_selectgroup(GROUP_ROTOR);
        for gi = 1:(2*p), mi_selectgroup(GROUP_ROTOR+gi); end
        mi_moverotate(0, 0, stepMechDeg);
        mi_clearselected;
    end
end

%% ---------------------- Back-EMF (circular central diff) --
dlamA_dth = (circshift(lamA,-1) - circshift(lamA,1)) / (2*dtheta_rad);
dlamB_dth = (circshift(lamB,-1) - circshift(lamB,1)) / (2*dtheta_rad);
dlamC_dth = (circshift(lamC,-1) - circshift(lamC,1)) / (2*dtheta_rad);
emfA = -dlamA_dth * omega_mech;
emfB = -dlamB_dth * omega_mech;
emfC = -dlamC_dth * omega_mech;

%% ---------------------- Save results ----------------------
results = struct;
results.meta = struct('p',p,'stepElecDeg',stepElecDeg,'stepMechDeg',stepMechDeg, ...
                      'rpm_mech',rpm_mech,'omega_mech',omega_mech, ...
                      'mode',mode,'Id',Id,'Iq',Iq,'theta_e_star_deg',theta_e_star_deg, ...
                      'rotationMechDeg',rotationMechDeg);
results.angle_mech = theta_mech(:);
results.angle_elec = (theta_e_star_deg + theta_elec(:));
results.torque     = T_mech;
results.lambda     = struct('A',lamA,'B',lamB,'C',lamC);
results.emf        = struct('A',emfA,'B',emfB,'C',emfC,'angle_mech',theta_mech(:));

save(outMat,'results');
fprintf('\nSimulation complete. Results saved to %s\n', outMat);

%% ---------------------- Cleanup temp files ----------------
closefemm;
keepTemp = false;                               % set true to keep tmp files
if ~keepTemp
    if exist(tmpFemFile,'file'), delete(tmpFemFile); end
    [tp,tn,~] = fileparts(tmpFemFile);
    fans = fullfile(tp,[tn '.ans']); if exist(fans,'file'), delete(fans); end
    fbak = fullfile(tp,[tn '.bak']); if exist(fbak,'file'), delete(fbak); end
end
