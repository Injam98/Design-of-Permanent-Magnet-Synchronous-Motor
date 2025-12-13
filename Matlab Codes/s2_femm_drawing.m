%% Run analytical model
run s1_geometry_calc.m

%% ========================= USER PARAMETERS ===============================
% FEMM installation path
femm_path = 'C:\femm42\mfiles';   % Adjust if FEMM is installed elsewhere

% FEMM simulation file name
femm_filename = 'PMSM_1.fem';

% FEMM problem settings
mesh_tolerance = 1e-8;            % Numerical tolerance
min_angle      = 30;              % Minimum element angle [deg]

% Mechanical shaft diameter (mm)
D_shaft = 55;                     % Rotor shaft diameter
R5 = D_shaft/2;                   % Shaft radius (currently not used, but defined)

% Slot copper layer radial offsets (mm)
offset_bottom = 10;               % Bottom (outer) copper label offset from x_slot(6)
offset_top    = 30;               % Top (inner) copper label offset from x_slot(6)

% Electrical angle for dq -> abc current assignment
theta_e = deg2rad(77);            % Electrical angle in radians (set as needed)

% Group IDs (used for geometry and post-processing)
GROUP_BOUNDARY = 1000;            % Outer boundary group
GROUP_SLOTS    = 1001;            % Stator slot group
GROUP_ROTOR    = 10;              % Rotor (yoke + magnets) group base
GROUP_STATOR   = 1;               % Stator iron group
COPPER_GROUP   = 2;               % Winding group (coils)

%% Generate femm script, Define Materials, Boundary and Circuit Property
addpath(femm_path);               % adding the path where the femm functions are located
currentFolder = fileparts(mfilename('fullpath'));
openfemm();                       % Open FEMM software
main_maximize;
newdocument(0);                   % Create a new FEMM document for a 2D magnetic problem
mi_probdef(0, 'millimeters', 'planar', mesh_tolerance, (l_prime*1000), min_angle, 0);
mi_smartmesh(0);
mi_addboundprop('outer', 0, 0, 0, 0, 0, 0, 0, 0, 0); % Boundary condition for 'Az'

%% ========================= Draw Stator ==================================
R_si = (D_si/2)*1000; % Inner radius of the stator
R_so = (D_so/2)*1000; % Outer radius of the stator

% Single slot dimensions
H1  = (h1+h2+h3+h4+h6)*1000; % Total height
H11 = h1*1000;               % Foot height
H13 = (h3+h5)*1000;          % Straight height
H14 = h4*1000*0.55;          % Coil separation height
B11 = b1*1000;               % Slot opening
B12 = b4*1000;               % Slot lower width
B13 = b5*1000;               % Slot higher width

phi_s = 2*pi/N_slot; % Dividing stator slot positions
x_slot(1) = R_si*cos(phi_s/2); y_slot(1) = -R_si*sin(phi_s/2);
x_slot(2) = sqrt(R_si^2-B11^2/4); y_slot(2) = -B11/2;
x_slot(3) = x_slot(2)+H11; y_slot(3) = y_slot(2);
x_slot(4) = x_slot(2)+H1-B13/2-H13; y_slot(4) = -B12/2;
x_slot(5) = x_slot(4)+H13; y_slot(5) = -B13/2;
x_slot(6) = x_slot(2)+H1; y_slot(6) = 0;
x_slot(7) = x_slot(5); y_slot(7) = -y_slot(5);
x_slot(8) = x_slot(4); y_slot(8) = -y_slot(4);
x_slot(9) = x_slot(3); y_slot(9) = -y_slot(3);
x_slot(10) = x_slot(2); y_slot(10) = -y_slot(2);
x_slot(11) = x_slot(1); y_slot(11) = -y_slot(1);
x_slot(12) = x_slot(4)+((H14*(x_slot(5)-x_slot(4)))/ ...
sqrt((x_slot(5)-x_slot(4))^2+(y_slot(5)-y_slot(4))^2));
y_slot(12) = y_slot(4)+((H14*(y_slot(5)-y_slot(4)))/ ...
sqrt((x_slot(5)-x_slot(4))^2+(y_slot(5)-y_slot(4))^2));
x_slot(13) = x_slot(12); y_slot(13) = -y_slot(12);
x_slot(14) = R_so*cos(phi_s/2); y_slot(14) = -R_so*sin(phi_s/2);
x_slot(15) = R_so*cos(phi_s/2); y_slot(15) = R_so*sin(phi_s/2);

% Add nodes
for n = 1: length(x_slot)
    mi_addnode(x_slot(n),y_slot(n));
end

% Connect arcs
mi_addarc(x_slot(1),y_slot(1),x_slot(2),y_slot(2),...
    rad2deg(phi_s/2-atan2(B11/2,sqrt(R_si^2-B11^2/4))),5);
mi_addarc(x_slot(10),y_slot(10),x_slot(11),y_slot(11),...
    rad2deg(phi_s/2-atan2(B11/2,sqrt(R_si^2-B11^2/4))),5);
mi_addarc(x_slot(5),y_slot(5),x_slot(6),y_slot(6),90,5);
mi_addarc(x_slot(6),y_slot(6),x_slot(7),y_slot(7),90,5);
mi_addarc(x_slot(14),y_slot(14),x_slot(15),y_slot(15),rad2deg(phi_s),5);

% Connect lines
for n = 2:4
    mi_addsegment(x_slot(n),y_slot(n),x_slot(n+1),y_slot(n+1));
end
for n = 9:-1:7
    mi_addsegment(x_slot(n+1),y_slot(n+1),x_slot(n),y_slot(n));
end
mi_addsegment(x_slot(4),y_slot(4),x_slot(8),y_slot(8));
mi_addsegment(x_slot(12),y_slot(12),x_slot(13),y_slot(13));

% Selecting every node of the slot
for n = 1:15
    mi_selectnode(x_slot(n),y_slot(n));
end
% Selecting every arc of the slot
mi_selectarcsegment(x_slot(1),y_slot(1));
mi_selectarcsegment(x_slot(10),y_slot(10));
mi_selectarcsegment(x_slot(5),y_slot(5));
mi_selectarcsegment(x_slot(7),y_slot(7));

% Selecting lines of one side of the slot
for n = 2:5
    mi_selectsegment((x_slot(n)+x_slot(n+1))/2,(y_slot(n)+y_slot(n+1))/2);
end

% Selecting lines of opposite side of the slot
for n = 9:-1:6
    mi_selectsegment((x_slot(n)+x_slot(n+1))/2,(y_slot(n)+y_slot(n+1))/2);
end

% Selecting segments
mi_selectsegment((x_slot(4)+x_slot(8))/2,(y_slot(4)+y_slot(8))/2);
mi_selectsegment((x_slot(12)+x_slot(13))/2,(y_slot(12)+y_slot(13))/2);
mi_setgroup(GROUP_SLOTS);
mi_clearselected;

mi_selectarcsegment(x_slot(15),y_slot(15));
mi_setarcsegmentprop(5, 'outer', 0, GROUP_BOUNDARY);
mi_selectnode(x_slot(14),y_slot(14));
mi_selectnode(x_slot(15),y_slot(15));
mi_setgroup(GROUP_BOUNDARY);
mi_clearselected;

%rotate = pi/N_slot;        % Rotate stator with windings 
%mi_selectgroup(1001);
%mi_moverotate(0, 0,rad2deg(rotate));
%mi_clearselected;

mi_selectgroup(GROUP_BOUNDARY);
mi_selectgroup(GROUP_SLOTS);
mi_copyrotate(0, 0, (360/N_slot), N_slot-1);
mi_clearselected;

mi_zoomnatural();

%% ============================= Draw rotor ===============================

phi_r = pi/p;                               % Dividing rotor slot positions
R1 = ((D_ro/2)-h_PM)*1000;                  % Rotor radius(inner) without PM
R2 = (D_ro/2)*1000;                         % Rotor radius(outer) with PM
R3 = ((D_ro/2)-h_PM-h_yr)*1000;             % Rotor radius without PM and Yoke
R4 = ((D_ro/2)-(h_PM/2))*1000;              % Magnet node radius
%R5 is defined in the user parameter block (shaft radius)
s1 = 2*R1*asin((w_PM*1000)/(2*R1));         % Inner radius arc length
s2 = s1*(1+((h_PM*1000)/R1));               % Outer radius arc length
ang = s1/R1;
x_rot(1) = R1*cos(-phi_r/2);      y_rot(1) = R1*sin(-phi_r/2);    % Inner arc 1st coordinate
x_rot(2) = R1*cos(phi_r/2);       y_rot(2) = R1*sin(phi_r/2);     % Inner arc 2nd coordinate
x_rot(3) = R2*cos(-ang/2);        y_rot(3) = R2*sin(-ang/2);      % Magnet outer arc 1st coordinate
x_rot(4) = R2*cos(ang/2);         y_rot(4) = R2*sin(ang/2);       % Magnet outer arc 2nd coordinate
x_rot(5) = R1*cos(-ang/2);        y_rot(5) = R1*sin(-ang/2);      % Magnet inner arc 1st coordinate
x_rot(6) = R1*cos(ang/2);         y_rot(6) = R1*sin(ang/2);       % Magnet inner arc 2nd coordinate
x_rot(7) = R3*cos(-phi_r/2);      y_rot(7) = R3*sin(-phi_r/2);    % Yoke 1st node
x_rot(8) = R3*cos(phi_r/2);       y_rot(8) = R3*sin(phi_r/2);     % Yoke 2nd node
%x_rot(9) = R5*cos(-phi_r/2);      y_rot(9) = R5*sin(-phi_r/2);    % Rotor shaft 1st coordinates
%x_rot(10) = R5*cos(phi_r/2);      y_rot(10) = R5*sin(phi_r/2);    % Rotor shaft 2nd coordinates

for n=1:length(x_rot)
    mi_addnode(x_rot(n),y_rot(n));  % Add nodes
end

mi_addsegment(x_rot(3),y_rot(3), x_rot(5),y_rot(5));                       % Connect nodes along magnet height
mi_addsegment(x_rot(4),y_rot(4), x_rot(6),y_rot(6));                       % Connect nodes along magnet height
mi_addarc(x_rot(1),y_rot(1), x_rot(5),y_rot(5), rad2deg((phi_r-ang)/2), 5);          % Connect inner arc nodes
mi_addarc(x_rot(2),y_rot(2), x_rot(6),y_rot(6), rad2deg((phi_r-ang)/2), 5);          % Connect inner arc nodes
mi_addarc(x_rot(3),y_rot(3), x_rot(4),y_rot(4), rad2deg(ang), 5);          % Connect outer arc nodes
mi_addarc(x_rot(5),y_rot(5), x_rot(6),y_rot(6), rad2deg(ang), 5);          % Connect outer arc nodes
mi_addarc(x_rot(7),y_rot(7), x_rot(8),y_rot(8), rad2deg(phi_r), 5);        % Connect yoke arc nodes
%mi_addarc(x_rot(9),y_rot(9), x_rot(10),y_rot(10), rad2deg(phi_r), 10);    % Connect shaft arc nodes

for n = 1:length(x_rot)
   mi_selectnode(x_rot(n),y_rot(n));    % Select nodes
end

% Known radii/angles for arcs:
% 1↔5 : radius R1, angles -phi_r/2 to -ang/2
% 2↔6 : radius R1, angles +phi_r/2 to +ang/2
% 3↔4 : radius R2, angles -ang/2   to +ang/2
% 5↔6 : radius R1, angles -ang/2   to +ang/2
% 7↔8 : radius R3, angles -phi_r/2 to +phi_r/2

% Helper inline to select an arc via its angular midpoint on radius R:
% (theta_mid in radians)
selArc = @(R,theta_mid) mi_selectarcsegment(R*cos(theta_mid), R*sin(theta_mid));

% Mid-angles
theta15 = ( -phi_r/2 + -ang/2 )/2;          % arc 1–5
theta26 = (  phi_r/2 +  ang/2 )/2;          % arc 2–6
theta34 = 0;                                % arc 3–4 (symmetric around x-axis)
theta56 = 0;                                % arc 5–6 (symmetric)
theta78 = 0;                                % arc 7–8 (symmetric)

% Select arc segments (points ON the arcs)
selArc(R1, theta15);                         % arc 1–5
selArc(R1, theta26);                         % arc 2–6
selArc(R2, theta34);                         % arc 3–4
selArc(R1, theta56);                         % arc 5–6
selArc(R3, theta78);                         % arc 7–8

% Select the straight magnet segments via chord midpoints:
mi_selectsegment((x_rot(3)+x_rot(5))/2, (y_rot(3)+y_rot(5))/2);
mi_selectsegment((x_rot(4)+x_rot(6))/2, (y_rot(4)+y_rot(6))/2);

% Group and copy
mi_setgroup(GROUP_ROTOR);
mi_selectgroup(GROUP_ROTOR);
mi_copyrotate(0, 0, 360/(2*p), (2*p-1));
mi_clearselected;

%{
mi_selectarcsegment(x_rot(8),y_rot(8));
mi_setgroup(14);
mi_selectgroup(14);
mi_copyrotate(0, 0, 360/(2*p), (2*p-1));
mi_clearselected;
%}

%% ===================== Define and add materials =========================

mi_getmaterial('Air');                                                      % Airgap material
mi_getmaterial('Copper');                                                   % Winding material
mi_addmaterial('NdFeB40', mu_rec, mu_rec, H_c, 0, 0, 0, 0, 1, 0, 0);        % Magnet material
mi_getmaterial('416 Stainless Steel');                                      % Shaft material
mi_addmaterial('M800-50A', 0, 0, 0, 0, 0, 0.5, 0, k_Fe, 0, 0, 0, 0, 0);     % Stator material

% Stator B-H values
BH_stator = [B_m
          H_m];
for nbh = 1:length(BH_stator)
    mi_addbhpoint('M800-50A', BH_stator(1,nbh),BH_stator(2,nbh));
end

%----------------------- Assign Airgap Materials -------------------------%
mi_addblocklabel((R2+((delta*1000)/2)), 0);
mi_selectlabel((R2+((delta*1000)/2)), 0);
mi_setblockprop('Air', 0, 0, 0, 0, 0, 0);
mi_clearselected;

%-------------------------Assign Stator Iron Material---------------------% 
mi_addblocklabel(R_si/2+H1/2+R_so/2, 0);
mi_selectlabel(R_si/2+H1/2+R_so/2, 0);
mi_setblockprop('M800-50A', 0, 0, 0, 0, GROUP_STATOR, 0);
mi_clearselected;

%-----------------------Assign Rotor Iron Materials-----------------------%
mi_addblocklabel((R1-((R1-R3)/2)), 0);  % Rotor yoke material node
mi_selectlabel((R1-((R1-R3)/2)), 0);
mi_setblockprop('M800-50A', 0, 0, 0, 0, GROUP_ROTOR, 0);
mi_clearselected;

mi_addblocklabel((R3/2), 0);            % Rest of the rotor material node
mi_selectlabel((R3/2), 0);
mi_setblockprop('Air', 0, 0, 0, 0, 0, 0);
%mi_setblockprop('M800-50A', 0, 0, 0, 0, 10, 0);
mi_clearselected;

%{
%-----------------------------Rotor Shaft---------------------------------% 
mi_addblocklabel(0,0);
mi_selectlabel(0,0);
mi_setblockprop('416 Stainless Steel', 0, 0, 0, 0, 9, 0);
mi_clearselected;
%}

%---------------------Assign Rotor Magnet Materials-----------------------%
% Place 2*p magnets, centers 360/(2p) apart, alternating out/in
x_mag = zeros(1, 2*p);
y_mag = zeros(1, 2*p);

for i = 1:(2*p)
    theta_deg = (i-1) * (360/(2*p));          % 0, 360/(2p), ... up to <360
    x_magn = R4 * cosd(theta_deg);
    y_magn = R4 * sind(theta_deg);

    % Alternate magnetization: outward for i=1, inward for i=2, etc.
    magndir_deg = theta_deg + 180 * mod(i-1, 2);   % add 180° on every other magnet

    mi_addblocklabel(x_magn, y_magn);
    mi_selectlabel(x_magn, y_magn);
    mi_setblockprop('NdFeB40', 0, 0, 0, magndir_deg, GROUP_ROTOR + i, 0);
    mi_clearselected;

    x_mag(i) = x_magn; 
    y_mag(i) = y_magn;
end


%%=============== Define Circuit and Assign Windings ======================

N_c = round(z_Qs/(2*a));                                            % Number of conductor in half slot

% Inverse Park transform (dq -> abc) using theta_e
Ia = I_d*cos(theta_e) - I_q*sin(theta_e);
Ib = I_d*cos(theta_e - 2*pi/3) - I_q*sin(theta_e - 2*pi/3);
Ic = I_d*cos(theta_e + 2*pi/3) - I_q*sin(theta_e + 2*pi/3);

Ia = sqrt(2)*Ia;
Ib = sqrt(2)*Ib;
Ic = sqrt(2)*Ic;

fprintf('ia = %.4f A\nib = %.4f A\nic = %.4f A\n', Ia, Ib, Ic);
mi_addcircprop('A', Ia, 1);  wind(1).name = 'A';                            % create Phase A circuit with current Ia
mi_addcircprop('B', Ib, 1);  wind(2).name = 'B';                            % create Phase B circuit with current Ib
mi_addcircprop('C', Ic, 1);  wind(3).name = 'C';                            % create Phase C circuit with current Ic

% -------- Single copper group id for ALL phases ------------------------------ 
% (COPPER_GROUP defined in user parameter block above)

% -------- Slot mapping (each entry = 'Top|Bottom') --------------------------- 
% Case encodes direction of that conductor: uppercase => +N_c, lowercase => -N_c.
mapStr = ["AA","bA","bb","Cb","CC","aC","aa","Ba","BB","cB","cc","Ac", ...
          "AA","bA","bb","Cb","CC","aC","aa","Ba","BB","cB","cc","Ac", ...
          "AA","bA","bb","Cb","CC","aC","aa","Ba","BB","cB","cc","Ac", ...
          "AA","bA","bb","Cb","CC","aC","aa","Ba","BB","cB","cc","Ac"];

% -------- Geometry helpers --------------------------------------------------- 
xslt0 = x_slot(6) - offset_bottom;   % radius for BOTTOM (outer) layer label
xslt1 = x_slot(6) - offset_top;      % radius for TOP (inner) layer label
phi_s  = 2*pi/N_slot;                % mechanical slot pitch [rad]

% -------- Helper: phase letter and signed turns from symbol ------------------
% Returns phaseName ('A'|'B'|'C') and turns (+/-N_c) based on case.
getPhaseAndTurns = @(sym) ...
    deal( upper(sym), ( (sym==upper(sym))*2 - 1 ) * N_c );                  % +N_c if uppercase, -N_c if lowercase

% -------- Place block labels for every slot (double-layer) -------------------
for n = 1:N_slot
    theta_slt = -(n*phi_s - (2*pi)/N_slot);   % your slot angular placement convention
    pair = char(mapStr(n));                   % two chars: [Top, Bottom], e.g., 'Aa', 'BB', etc.

    % ----- TOP conductor (INNER layer) ---------------------------------------
    topSym = pair(1);                                                       % symbol for top conductor in this slot
    % sanity check (helps catch typos in mapStr)
    if ~ismember(topSym, ['A','B','C','a','b','c'])
        error(['mapStr(%d) top symbol ''%s'' is invalid. ' ...
            'Use only A/B/C/a/b/c.'], n, topSym);
    end
    [topPhase, topTurns] = getPhaseAndTurns(topSym);                        % resolve phase name and +/-N_c
    xt = cos(theta_slt)*xslt1;  yt = sin(theta_slt)*xslt1;                  % coordinates for TOP label (inner)
    mi_addblocklabel(xt, yt);                                               % add label at TOP position
    mi_selectlabel(xt, yt);                                                 % select the just-added label
    mi_setblockprop('Copper', 0, 0, topPhase, 0, COPPER_GROUP, topTurns);   % Block property
    mi_clearselected;                                                       % clear selection

    % ----- BOTTOM conductor (OUTER layer) ------------------------------------
    botSym = pair(2);                                  % symbol for bottom conductor in this slot
    if ~ismember(botSym, ['A','B','C','a','b','c'])
        error(['mapStr(%d) bottom symbol ''%s'' is invalid. ' ...
            'Use only A/B/C/a/b/c.'], n, botSym);
    end
    [botPhase, botTurns] = getPhaseAndTurns(botSym);                        % resolve phase name and +/-N_c
    xb = cos(theta_slt)*xslt0;  yb = sin(theta_slt)*xslt0;                  % coordinates for BOTTOM label (outer)
    mi_addblocklabel(xb, yb);                                               % add label at BOTTOM position
    mi_selectlabel(xb, yb);                                                 % select the just-added label
    mi_setblockprop('Copper', 0, 0, botPhase, 0, COPPER_GROUP, botTurns);   % signed turns for BOTTOM conductor
    mi_clearselected;                                                       % clear selection

    % ----- Debug printout (shows phase & signed turns for each layer) --------
    fprintf('Slot %2d: TOP=%s (%+dT)  BOTTOM=%s (%+dT)  [Group=%d]\n', ...
        n, topPhase, topTurns, botPhase, botTurns, COPPER_GROUP);
end


%% Save the geometry and analyze
mi_saveas(femm_filename);
mi_analyze;
mi_loadsolution;
mo_groupselectblock(GROUP_ROTOR);
for i=1:(2*p)
    mo_groupselectblock(GROUP_ROTOR + i);
end
T_FEMM = mo_blockintegral(22);
fprintf('Torque (WST) = %.3f N·m\n', T_FEMM);
mo_showdensityplot(1,0,0,2,'bmag');
%mo_close;
%{
propsA = mo_getcircuitproperties('A');  IA = propsA(1);  VA = propsA(2);  PA = propsA(3);
propsB = mo_getcircuitproperties('B');  IB = propsB(1);  VB = propsB(2);  PB = propsB(3);
propsC = mo_getcircuitproperties('C');  IC = propsC(1);  VC = propsC(2);  PC = propsC(3);

fprintf('Solved currents: A=%.2f  B=%.2f  C=%.2f A\n', IA, IB, IC);

%mi_addmaterial(’matname’, mu x, mu y, H c, J, Cduct, Lam d, Phi hmax, lam fill, LamType, Phi hx, Phi hy, nstr, dwire)
%mi_setblockprop(’blockname’, automesh, meshsize, ’incircuit’, magdir, group, turns)
%mi_setsegmentprop(’propname’, elementsize, automesh, hide, group)
%}
