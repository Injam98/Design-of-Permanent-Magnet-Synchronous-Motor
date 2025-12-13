%% torque_vs_thetae.m : rotor fixed, sweep current injection angle
% Requires:
% - s1_geometry_calc.m (analytical design) already run
% - s2_FEMM_Geometry.m (FEMM model build) already run
% so that: femm_filename, p, I_d, I_q, GROUP_ROTOR exist in the workspace.

%% USER PARAMETERS (edit if needed)
stepElecDeg      = 1;                           % electrical angle step [deg]
rotationDeg      = 164;                         % electrical span to scan [deg]

mat_output_file  = 'torque_vs_thetae_peak.mat';  % .mat file to save results
temp_fem_file    = 'temp.fem';                  % temporary FEMM model file
temp_ans_file    = 'temp.ans';                  % temporary FEMM solution file

%% FEMM setup and open model
openfemm();
main_maximize;
opendocument(femm_filename);                    % FEMM file from Script 2
mi_smartmesh(0);                                % coarse mesh to be quick
mi_saveas(temp_fem_file);                      % work on a temp copy

theta_e_deg_vec = -15:stepElecDeg:(rotationDeg-stepElecDeg);  % sweep range
nSteps = numel(theta_e_deg_vec);               % number of samples
Tq = zeros(1,nSteps);                          % torque array

% (Circuits exist: 'A','B','C')

for k = 1:nSteps                               % loop over theta_e
    theta_e = deg2rad(theta_e_deg_vec(k));     % current injection angle [rad]
    
    % dq -> abc using I_d, I_q from Script 1
    Ia = I_d*cos(theta_e) - I_q*sin(theta_e);                     % phase a
    Ib = I_d*cos(theta_e - 2*pi/3) - I_q*sin(theta_e - 2*pi/3);   % phase b
    Ic = I_d*cos(theta_e + 2*pi/3) - I_q*sin(theta_e + 2*pi/3);   % phase c
    
    Ia = sqrt(2) * Ia;
    Ib = sqrt(2) * Ib;
    Ic = sqrt(2) * Ic;

    mi_modifycircprop('A', 1, Ia);             % set phase-A current
    mi_modifycircprop('B', 1, Ib);             % set phase-B current
    mi_modifycircprop('C', 1, Ic);             % set phase-C current
    
    mi_createmesh();
    mi_analyze();
    mi_loadsolution();
    mo_clearblock;                             % clear selection
    
    % Select rotor yoke and magnets (GROUP_ROTOR from Script 2)
    mo_groupselectblock(GROUP_ROTOR);          % rotor yoke
    for m = 1:(2*p)
        mo_groupselectblock(GROUP_ROTOR + m);  % magnets
    end
    
    Tq(k) = mo_blockintegral(22);              % torque [N·m]
    fprintf('theta_e=%6.1f°  T=%.3f N·m\n', theta_e_deg_vec(k), Tq(k)); % log
end

[Tmax, idxMax] = max(Tq);                      % peak torque and index
theta_e_star_deg = theta_e_deg_vec(idxMax);    % best injection angle [deg]
fprintf('Max Torque = %.3f N·m at theta_e* = %.1f° (rotor fixed)\n', Tmax, theta_e_star_deg);

save(mat_output_file, ...
     'theta_e_deg_vec','Tq','theta_e_star_deg','Tmax');  % save

figure;     % quick plot
plot(theta_e_deg_vec, Tq,'LineWidth',1.5); grid on;
xlabel('\theta_e (deg)'); ylabel('Torque (N·m)');
title('Torque vs \theta_e');

delete(temp_fem_file);
delete(temp_ans_file);                          % cleanup
