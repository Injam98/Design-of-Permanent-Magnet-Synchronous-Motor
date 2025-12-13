clc; clear; close all

%% 1. Initial Data
P = 160000;                 % Shaft power in watt
rpm_mech = 2.5000e3;        % Speed in RPM
n = rpm_mech/60;            % Speed in rps
T = P/(2*pi*n);             % Torque in Nm
U = 690;                    % Line-to-line voltage in volt
U_sph = U/sqrt(3);          % Stator phase voltage in V
m = 3;                      % Number of phases
p = 4;                      % Number of pole pairs
f = n*p;                    % Frequency in Hz
omega = 2*pi*f;             % Angular frequency in 1/s
eta = 0.95;                 % Efficiency, maximum value must be 1
pf = 0.91;                  % value of cos(phi)
H_c = 800000;               % Coercivity of PM in A/m
%H_c = 883310;
mu_0 = 4*pi*10^-7;          % Permeability of vacuum in V.s/A.m
B_r = 1.05;                 % Remanence flux density of permanent magnet in T
mu_rec = B_r/(mu_0*H_c);    % Permeability of the permanet magnet material
theta = 80;                 % Temperature rise in the machine windings in K
sigma_Cu20C = 57e6;         % Conductivity of copper at 20 degrees C in S/m
alpha_Cu = 3.81e-3;         % Temperature coefficient of resistivity for copper in 1/K
k_Fe = 0.97;                % Space factor of the stator core
rho_Fe = 7600;              % Density of iron in kg/m^3
rho_PM = 7500;              % Density of the permanent magnet material in kg/m^3
rho_Cu = 8960;              % Density of copper in kg/m^3
P_15 = 6.6;                 % Specific loss of lamination material at 1.5 T and 50 Hz in W/kg
B_m = (0:0.1:2);            % Magnetic flux density of lamination material in T
H_m = [0, 84.5, 107 121 133 145 156 168 180 194 209 228 254 304 402 660 ...
    1480 3710 7300 15000 30000];    % Magnetic field intensity of lamination material in A/m
B_y = (0:0.1:1.9)';                 % Yoke flux density
c = [0.72 0.72 0.72 0.72 0.72 0.72 0.71 0.70 0.67 0.63 0.57 0.48 0.40 0.33...
    0.26 0.20 0.17 0.16 0.15 0.14]';% Coefficient

% Cooling channels (no cooling ducts)
n_v = 0; b_v = 0; b_ve = 0;         % No cooling channels

% Winding / slotting design
q = 2;                      % Number of slots per pole per phase (defined)
W_taup = 5/6;               % Winding pitch

% Air-gap flux design
B_1peak = 0.95;             % Fundamental air gap flux peak density in T (chosen)
alpha_PM = 0.8;             % Effective relative magnet width (chosen)

% Tooth flux density design
B_dapp = 1.6;               % Apparent flux density of stator tooth (chosen)

% Current density & slot fill
J_s = 4.5e6;                % Stator current density in A/m^2 chosen from table 6.2
k_Cus = 0.63;               % Copper filling factor of the slot

% Slot geometry (chosen)
b1 = 0.003;                 % Slot opening width
h1 = 0.001;                 % Slot height section 1
h2 = 0.002;                 % Slot height section 2
h3 = 0.005;                 % Slot height section 3
h6 = 0.0005;                % Slot lip height
h_prime = 0.0005;           % Additional slot height used in leakage calc

% h5 iteration settings (slot copper area)
h5 = 0.001;                 % Initial guess for h5
tolerance_h5 = 1e-9;        % Convergence criteria for h5
max_iterations_h5 = 1000;   % Safety limit for h5 loop

% Yoke flux densities
B_ys = 1.3;                 % Stator yoke flux density (from Table 6.1)
B_yr = 1.3;                 % Rotor yoke flux density (from Table 6.1)

% End winding leakage parameters
l_ew = 0.025;               % Chosen end-winding overhang length
lambda_lew = 0.518;         % Chosen from table 4.2
lambda_W = 0.138;           % Leakage coefficient

% Mechanical / PM loss parameters
k_rho = 10;                 % Mechanical loss coefficient (chosen)
sigma_PM = 670000;          % PM conductivity varies with material property

% Power matching (delta_l) iteration parameters
tolerance_P = 1e-5;         % Tolerance for power matching
step_size = 1e-5;           % Step size for delta_l iteration

%% 2. Tangential Stress
sigma_ftan = 33500*pf;      % Average Tangential stress in Pa from table 6.3

%% 3. Rotor Size
V_r = T/(2*sigma_ftan);                     % Rotor volume in m^3 from equation 6.2
chi = (pi*p^0.5)/(4*p);                     % Ratio of equivalent core length and air-gap diameter from table 6.5
D_ro = round(((4*V_r)/(pi*chi))^(1/3),2);   % Rotor outer diameter in m
l_prime = chi*D_ro;                         % Core length in m

%% 4. Air gap and Core length
delta = (0.18+(0.006*P^0.4))/1000;  % Air gap length in m from surface of PM to stator inner surface
delta = delta+1e-3;                 % Added surface magnet supporting band of 1mm
D_si = D_ro+(2*delta);              % Inner diameter of stator in m
l = l_prime-(n_v*b_ve)-(2*delta);   % "l" = Core length in m

%% 5. Stator Winding
N_slot = 2*p*m*q;           % Number of stator slots
tau_u = (D_si*pi)/N_slot;   % Stator slot pitch
tau_p = (D_si*pi)/(2*p);    % Stator pole pitch

%% 6. Air-gap flux density and linear current density
B_max = ((pi*B_1peak)/(4*sin((alpha_PM*pi)/2)));    % Maximum value of rectangular flux density
A = (sigma_ftan/(0.5*B_1peak*pf*sqrt(2)));          % RMS linear current density in A/m, reference value is in table 6.3 for reference

%% 7. Number of coil turns in a phase winding
E_pm = U/sqrt(3);                                       % Permanent magnet induced voltage
k_w1 = (2*(sin((pi*W_taup)/2))*(sin(pi/(2*m))))/...
    ((N_slot/(m*p))*(sin((pi*p)/N_slot)));              % Winding factor for fundamental
N_turn = (sqrt(2)*E_pm)/...
    (omega*k_w1*alpha_PM*B_max*tau_p*l_prime);          % Number of coil turns in series in a phase winding

%% 8. Number of conductors in a slot
a = 1;                          % Number of parallel branches (Chosen)
% Compute z_Q and check if it's even
while true
    z_Q = (2*a*m*N_turn)/N_slot;          % Calculated number of conductor in a slot
    if mod(floor(z_Q), 2) == 0  % Check if the rounded value is even
        z_Qs = round(z_Q);      % Set z_Q to the nearest even integer (Selected number of conductor in a slot)
        break;                  % Exit loop if z_Q is even
    else
        a = a + 1;              % Increase 'a' if z_Q is not even
    end
end

N_turn = (N_slot*z_Qs)/(2*a*m);           % Compute N (number of coil turns in the phase winding)

%% 9. New Maximum air-gap flux density
B_maxnew = (z_Q/z_Qs)*B_max;    % New air-gap flux density for selected z_Qs value in T
B_max = B_maxnew;

%% 10. Width of the stator slot
b_d = (l_prime*tau_u*B_max)/(k_Fe*B_dapp*(l-(n_v*b_v)));    % Tooth width in m

%% 11. Stator slot dimensions
I_s = P/(m*eta*U_sph*pf);   % Estimated stator current in A
S_cs = I_s/(a*J_s);         % Area of the conductor in a stotor slot
S_Cus = (z_Qs*S_cs)/k_Cus;  % Wound area of the stator slot
b4 = ((pi*(D_si+(2*(h1+h2))))/N_slot)-b_d;
b_4c = b4+((2*pi*h3)/N_slot)-(2*h6);

% Initial guess for h_5
S_Cus_target = S_Cus;
iteration = 0;

while iteration < max_iterations_h5
    b_5c = b_4c+(2*pi*h5)/N_slot;                         % Calculate b_5c using the first equation
    S_Cus = ((b_4c+b_5c)/2)*h5+(pi/8)*b_5c^2;             % Calculate S_Cus using the second equation
    
    % Check for convergence
    error = abs(S_Cus - S_Cus_target);
    if error < tolerance_h5
        break;  % Solution found
    end
    h5 = h5 * (S_Cus_target / S_Cus);                     % Adjust h_5 using a small step (simple numerical adjustment)
    iteration = iteration + 1;
end
b5 = b_5c+(2*h6);
h4 = h5+(b_5c/2);
S_slot = (b1*h1)+(h2*((b4/2)+(b1/2)))+(h3*(b4+...
    ((pi*h3)/N_slot)))+(h5*((b4+b5)/2))+((pi/8)*b5^2);  % Total area of the slot in m^2
h_d = h1+h2+h3+h4+h6;                                   % Slot depth of stator

%% 12. Magnetic voltage over the tooth
% Plot the B-H curve
figure();
plot(H_m,B_m, 'Linewidth', 1); grid on;
xlabel('Magnetic Field Intensity, H (A/m)'); ylabel('Flux Density, B (T)');
title('B-H curve of Iron (M800-50)');

H_interp = interp1(B_m, H_m, B_dapp, 'spline');                     % Find H for corresponding B
B_d = B_dapp-((((l_prime*tau_u)/(k_Fe*l*b_d))-1)*mu_0*H_interp);    % Find B for tooth
H_d = interp1(B_m, H_m, B_d, 'spline');                             % Find H for tooth from corresponding B

U_mds = H_d*(h3+h5);      % Magnetic voltage of tooth tip and rounded part ommitted due to small

%% 13. Checking the saturation factor
% In the case of a permanent magnet motor, the saturation of the teeth does
% not influence the waveform of the air-gap flux density, in other words, 
% α_PM remains constant

%% 14. Height of stator and rotor yokes and permanent magnets. Magnetic voltage of air gap, stator and rotor yokes, and permanent magnet
% Magnetic voltage of air gap
kappa = (2/pi)*((atan(b1/(2*delta)))-(((2*delta)/b1)*log(sqrt(1+(b1/(2*delta))^2))));
k_c1 = tau_u/(tau_u-(kappa*b1));   % Carter factor
delta_e = k_c1*delta;              % Equivalent air gap
U_mdeltae = (B_max*delta_e)/mu_0;  % Air gap magnetic voltage

% Stator and rotor yoke magnetic voltage
phi_m = alpha_PM*B_max*tau_p*l_prime;       % Air gap flux
h_ys = phi_m/(2*k_Fe*B_ys*(l-(n_v*b_v)));   % Stator yoke height
%h_ys = 0.1;
h_yr = phi_m/(2*k_Fe*B_yr*(l-(n_v*b_v)));   % Rotor yoke height

% Plot the B-c curve
figure();
plot(B_y,c, 'Linewidth', 1); grid on;
xlabel('Yoke Flux Density, B_y (T)'); ylabel('Coefficient, c');
title('Influence of maximum flux density of stator/rotor yoke on the definition of coefficient c, Fig-3.17');

H_ysmax = interp1(B_m, H_m, B_ys, 'spline');     % Find H for corresponding B
H_yrmax = interp1(B_m, H_m, B_yr, 'spline');     % Find H for corresponding B

c_s = interp1(B_y, c, B_ys, 'spline');     % Find correction factor c for corresponding B_ys;
c_r = interp1(B_y, c, B_yr, 'spline');     % Find correction factor c for corresponding B_yr;

D_ys = D_si+(2*(h1+h2+h3+h4+h6))+(2*h_ys);  % Stator yoke diameter
tau_ys = (pi*D_ys)/(2*p);                   % Pole pith on outer stator diameter
U_mys = c_s*H_ysmax*tau_ys;                 % Magnetic voltage of the stator yoke

% Permanent magnet magnetic voltage
B_PM = B_max;
h_PM = (U_mdeltae+U_mds+(U_mys/2)+((pi*c_r*H_yrmax*(D_ro-h_yr))/(4*p)))/...
    (H_c-((H_c*B_PM)/B_r)+((pi*c_r*H_yrmax)/(2*p)));                        % Permanent magnet height
w_PM = tau_p*alpha_PM;                                                      % Magnet width
U_mPM = (H_c/B_r)*B_PM*h_PM;

% The magnetic voltage of the rotor yoke
D_yr = D_ro-(2*h_PM);           % Average rotor yoke diameter
tau_yr = (pi*D_yr)/(2*p);       % Pole pitch on rotor yoke diameter
U_myr = c_r*H_yrmax*tau_yr;     % Rotor yoke magnetic voltage

%% 15. Outer stator and inner rotor diameter
D_so = D_ys;            % Stator outer diameter
D_ri = D_yr-(2*h_yr);   % Rotor inner diameter

%% 16. Total Magnetic voltage of the magnetic circuit
U_mtot = U_mdeltae+U_mds+U_mPM+(U_mys/2)+(U_myr/2);

%% 17. Stator Resistance
l_avg = (2*l)+(2.4*W_taup*tau_p)+0.1;           % Average length of a coil turn
sigma_Cu = sigma_Cu20C/(1+(theta*alpha_Cu));    % Conductivity of copper wire in 100 degrees
R=(N_turn*l_avg)/(sigma_Cu*a*S_cs);             % DC resistance of a phase winding

%% 18. Magnetizing Inductance
delta_eff = (U_mtot*delta_e)/U_mdeltae;                           % Effective air gap
L_md = (m/2)*(2/pi)*mu_0*l_prime*(1/(2*p))*(4/pi)*...
    (tau_p/delta_eff)*(k_w1*N_turn)^2;                            % d-axis inductance
L_mq = L_md;                                                       % Because of symmetrical

%% 19. Air-gap leakage inductance and reactance
alpha_u = (p*2*pi)/N_slot;                     % Slot angle

% First part of air gap leakage
k_delta1 = 0;
for k = 1:300
    numerator = sin((1+2*k*m)*W_taup*(pi/2))*(sin((1+2*k*m)*...
    q*(alpha_u/2))/(q*sin((1+2*k*m)*(alpha_u/2))));
    denominator = (1+2*k*m)*k_w1;
    value = (numerator/denominator)^2;
    k_delta1 = k_delta1+value; % Accumulate sum
end
clear numerator denominator value
%fprintf('The total sum is: %.20f\n', k_delta1);

% Second part of air gap leakage
k_delta2 = 0;
for k = -300:-1
    numerator = sin((1+2*k*m)*W_taup*(pi/2))*(sin((1+2*k*m)* ...
    q*(alpha_u/2))/(q*sin((1+2*k*m)*(alpha_u/2))));
    
    denominator = (1+2*k*m)*k_w1;
    
    value = (numerator/denominator)^2;
    k_delta2 = k_delta2+value; % Accumulate sum
end
clear numerator denominator value
%fprintf('The total sum is: %.20f\n', k_delta2);

sigma_deltas = k_delta1+k_delta2;     % Air gap leakage factor
L_delta = sigma_deltas*L_md;          % Air gap leakage inductance
X_delta = L_delta*2*pi*f;             % Air gap leakage reactance

%% 20. Slot leakage inductance and reactance
epsilon = 1-W_taup;
k_1 = 1-((9*epsilon)/16);
k_2 = 1-((3*epsilon)/4);

lambda_u = ((k_1*(h4-h_prime))/(3*b4))+(k_2*((h3/b4)+(h1/b1)+...
           (h2/(b4-b1))*log(b4/b1)))+(h_prime/(4*b4));
L_u = (4*m*mu_0*l_prime*N_turn^2*lambda_u)/N_slot;  % Slot leakage inductance
X_u = 2*pi*f*L_u;                                   % Slot leakage reactance

%% 21. Tooth tip leakage inductance and reactance
lambda_d = k_2*(((5*delta)/b1)/(5+((4*delta)/b1)));
L_sigmad = ((4*m)/N_slot)*mu_0*l_prime*lambda_d*N_turn^2;  % Tooth tip leakage inductance
X_sigmad = 2*pi*f*L_sigmad;                                % Tooth tip leakage reactance

%% 22. End winding leakage inductance and reactance
l_w = (l_avg/2)-l;
W_ew = l_w-(2*l_ew);
lambda_w = ((2*l_ew*lambda_lew)+(W_ew*lambda_W))/l_w;
L_w = ((4*m) / N_slot)*q*mu_0*l_w*lambda_w*N_turn^2;    % End winding Leakage inductance
X_w = 2*pi*f*L_w;                                       % End winding leakage reactance

%% 23. Synchronous inductance and reactance

L_ssigma = L_delta+L_u+L_sigmad+L_w;    % Stator leakage inductance
X_ssigma = 2*pi*f*L_ssigma;             % Stator leakage reactance
L_d = L_md+L_ssigma;                    % Direct axis synchronous inductance and reactance
X_d = 2*pi*f*L_d;                       % Direct axis synchronous reactance
L_q = L_d;                              % due to no saliency
X_q = X_d;                              % due to no saliency

%% 24. Losses (except stator losses)
V_primes = (pi/4)*(D_so^2-D_si^2)*l;            % Total volume of stator
V_ys = pi*l*((D_so/2)^2-((D_so/2)-h_ys)^2);     % Stator yoke volume
m_ys = V_ys*k_Fe*rho_Fe;                        % Stator yoke mass
V_slots = N_slot*S_slot*l;                      % Total volume of stator slot
V_ds = V_primes-V_ys-V_slots;                   % Teeth volume

% Total mass of teeth
m_ds = V_ds*k_Fe*rho_Fe;                % Only h_5 taken intor account for loss calculation due to low flux density of other parts, it is ignored
m_d = k_Fe*rho_Fe*N_slot*b_d*h5*l;

% Correction coefficients for core loss calculations Table 3.2
k_Fed = 2;
k_Fey = 1.5;

% Core loss in the stator yoke from eq 3.69 and 3.77
P_Feys = k_Fey*P_15*((B_ys/1.5)^2*m_ys*(f/50)^(3/2));

% Core loss in the tooth area
P_Feds = k_Fed*P_15*((B_d/1.5)^2*m_d*(f/50)^(3/2));

% Total iron loss
P_Fe = P_Feys+P_Feds;

% Mechanical losses (Windage and Ventilation) from experimental equation Eq 9.19 and table 9.2
v_r = pi*n*D_ro;
P_rho = k_rho*D_ro*v_r^2*(l+(0.6*tau_p));

% Effect of slot openings b_1 causes frequency f_PM on the rotor surface
f_PM = n*N_slot;

% Fictitious air gap for loss calculation
delta_PMEC = delta+(h_PM/(2*mu_rec));
b_e = kappa*b1;
k_CPM = tau_u/(tau_u-b_e);
u = (b1/(2*delta_PMEC))+sqrt(1+(b1/(2*delta_PMEC))^2);
beta = (1+u^2-(2*u))/(2*(1+u^2));
B_0 = beta*B_max;
k_v = sqrt(2*pi*f_PM*mu_rec*mu_0*(sigma_PM/2));
beta_v = (f_PM*2*pi)/(pi*D_ro*n);
a_Rv = (1/sqrt(2))*sqrt((sqrt(4+(beta_v/k_v)^4)+(beta_v/k_v)^2));
P_PMEC = (a_Rv/2)*(1+(tau_u/(2*l)))*(B_0/(mu_0*mu_rec))^2 * ...
    (k_v/sigma_PM)*pi*D_ro*alpha_PM*l*((k_v*sqrt(2))^2/beta_v^2);

% Additional losses
P_ex = 0.005*P;

%% 25. Rated load, stator current, stator resistive losses and total losses
while true
    % Stator resistive loss
    P_Cu = 3*I_s^2*R;
    
    % Sum of losses
    P_loss = P_Fe+P_Cu+P_rho+P_PMEC+P_ex;

    % Input power
    P_in = P+P_loss;

    % Initial values
    delta_l = 0;               % Initial guess for delta_l in radians
    P_in_given = P_in;

    % Iteration to find delta_l
    while true
        % Calculate P_in based on current delta_l from Eq. 7.214
        P_in = 3*(((U_sph*E_pm*sin(delta_l))/X_d)+ ...
               ((U_sph^2*(X_d-X_q)*sin(2*delta_l))/(2*X_d*X_q)));

        % Check if P_in matches the given value within tolerance
        if (P_in_given-P_in)<tolerance_P
            diff = (P_in_given-P_in);
            break;
        end

        % Adjust delta_l based on power difference
        if P_in<P_in_given
            delta_l = delta_l+step_size; % Increase delta_l if power is low
        else
            delta_l = delta_l-step_size; % Decrease delta_l if power is high
        end

        %diff = (P_in_given - P_in);
        %fprintf('Difference: %.5f W, delta: %.5f\n', diff, delta_l);
    end
    
    % Convert final delta_l to degrees
    delta_l_deg = rad2deg(delta_l);

    % Display results
    %fprintf(['Estimated Power: %.5f and Calculated Power: %.5f W at ' ...
    %   'Final delta_l: %.5f radians (%.5f degrees)\n Difference: %.5f'], ...
    %  P_in_given, P_in, delta_l, delta_l_deg, diff);
    
    % Direct-axis component of the stator current
    I_d = (((U/sqrt(3))*((X_q*cos(delta_l))-(R*sin(delta_l))))-(E_pm*X_q))/((X_d*X_q)+R^2);
    
    % Quadrature-axis component of the stator current
    I_q = (((U/sqrt(3))*((R*cos(delta_l))+(X_d*sin(delta_l))))-(E_pm*R))/((X_d*X_q) + R^2);
    
    % Stator current
    I_s_new = sqrt(I_d^2 + I_q^2);
    
    %fprintf('Loop\n');
    
    % Check if the difference between I_s and I_s_new is within 1%
    if abs((I_s_new-I_s)/I_s)<=0.01
        diffp = 100*abs((I_s_new-I_s)/I_s);
        I_s=I_s_new;
        fprintf('Difference: %.5f percent and final I_s is %.5f\n', diffp, I_s);
        break;
    else
        % Update I_s for the next iteration
        I_s = I_s_new;
    end
end
fprintf('I_s final: %.5f amp\n', I_s);

%% 26. Efficiency and power factor
eta_c = (P/P_in)*100;           % Effifiency of the machine
pf_c = P_in/(sqrt(3)*U*I_s);    % pf


%% 27. Masses of the active materials
D_ryi = D_ro-(2*h_PM);                             % Rotor core outer diameter
m_PM = ((D_ryi+D_ri)/2)*pi*alpha_PM*l*h_PM*rho_PM;  % PM mass
m_Cu = rho_Cu*(l+(2*l_w))*N_slot*z_Q*S_cs;               % Cu mass
m_yr = k_Fe*pi*l*rho_Fe*((D_ryi^2-D_ri^2)/4);     % Rotor yoke mass
m_tot = m_yr+m_ys+m_ds+m_Cu+m_PM;

%% 28. Load angle equation graph
i = 0:1:100;
delta_load = i*(pi/100);
T_delta_load = ((3*p)/(2*pi*f))*(((E_pm*U_sph*sin(delta_load))/X_d) +...
                ((U_sph^2/2)*sin(2*delta_load)*(X_q^-1 - X_d^-1)));

figure()
plot(rad2deg(delta_load), T_delta_load);
xlabel('Load Angle (deg)'); ylabel('Torque (N-m)');
title('Load characteristics without losses');
grid on;

%% 29. Summary of the machine design data
% Main Data
fprintf('\nMain Data\n');
fprintf('Input power, P_in = %.5f\n', P_in);
fprintf('Output power, P = %.5f\n', P);
fprintf('Efficiency, eta = %.2f\n', eta_c);
fprintf('PF = %.2f\n', pf_c);
fprintf('Input voltage, U = %.2f V\n', U);
fprintf('Input current, I_s = %.4f A\n', I_s_new);
fprintf('Input current, I_d = %.4f A\n', I_d);
fprintf('Input current, I_q = %.4f A\n', I_q);
fprintf('Load angle, delta_l = %.4f A\n', delta_l_deg);
fprintf('Number of turns, N = %d\n', N_turn);
fprintf('Carter factor, K_C1 = %.4f\n', k_c1);
fprintf('Magnet height, h_PM = %.4f m\n', h_PM);
fprintf('Slot number, Q = %d\n', N_slot);
fprintf('Number of slot per pole per phase, q = %d\n', q);
fprintf('Winding pitch, W_taup = %.4f\n', W_taup);
fprintf('Maximum flux density, B_max = %.4f\n', B_max);
fprintf('Machine length, l_prime = %.4f\n', l_prime);
fprintf('Stator inner diameter, D_s = %.4f\n', D_si);
fprintf('Stator outer diameter, D_se = %.4f m\n', D_so);
fprintf('Stator pole pitch, tau_p = %.4f\n', tau_p);
fprintf('Stator slot pitch, tau_u = %.4f\n', tau_u);
fprintf('Effective relative magnet width, alpha_PM = %.4f\n', alpha_PM);

fprintf('\nLosses\n');
fprintf('Copper loss, P_Cu = %.4f W\n', P_Cu);
fprintf('Iron Loss, P_Fe = %.4f W\n', P_Fe);
fprintf('Additional loss (0.5 * P), P_ex = %.4f W\n', P_ex);
fprintf('Losses including windage and ventilator, P_rho = %.4f W\n', P_rho);
fprintf('Mechanical Losses, P_MEC = %.4f W\n', P_PMEC);

fprintf('\nMagnetic Voltages\n');
U_m = (2*U_mPM) + (2*U_mds) + (2*U_mdeltae) + U_mys + U_myr;
fprintf('MMF, U_m = %.4f 1/A\n', U_m);
fprintf('PM magnetic voltage, U_mPM = %.4f 1/A\n', U_mPM);
fprintf('Stator yoke magnetic voltage, U_mys = %.4f 1/A\n', U_mys);
fprintf('Stator tooth magnetic voltage, U_mds = %.4f 1/A\n', U_mds);
fprintf('PM Rotor yoke magnetic voltage, U_myr = %.4f 1/A\n', U_myr);
fprintf('Air gap magnetic voltage, U_mdeltae = %.4f 1/A\n', U_mdeltae);
theta_PM = h_PM * H_c;
fprintf('theta_PM = %.4f 1/A\n', theta_PM);
ttheta_PM = 2*theta_PM;
fprintf('2*theta_PM = %.4f 1/A\n', ttheta_PM);

fprintf('\nSlot Dimensions\n');
fprintf('b_1 = %.4f m\n', b1);
fprintf('b_4 = %.4f m\n', b4);
fprintf('b_5 = %.4f m\n', b5);
fprintf('h_1 = %.4f m\n', h1);
fprintf('h_2 = %.4f m\n', h2);
fprintf('h_3 = %.4f m\n', h3);
fprintf('h_4 = %.4f m\n', h4);
fprintf('h_5 = %.4f m\n', h5);
fprintf('h_6 = %.4f m\n', h6);

fprintf('\nResistances and Reactances\n');
Z_N = U_sph / I_s;
fprintf('Z_N = %.4f Ohm', Z_N);
L_N = Z_N / omega;
fprintf('\nL_N = %.4f H\n', L_N);
fprintf('R = %.4f Ohm\n', R);
fprintf('X_d = %.4f H\n', X_d);
fprintf('L_d = %.4f H\n', L_d);
fprintf('L_md = %.4f H\n', L_md);
fprintf('L_q = %.4f H\n', L_q);
fprintf('L_mq = %.4f H\n', L_mq);
fprintf('X_ssigma = %.4f H\n', X_ssigma);
fprintf('L_ssigma = %.4f H\n', L_ssigma);
R_pu = R / Z_N;
fprintf('R_pu = %.4f Ohm\n', R_pu);
L_dpu = L_d / L_N;
fprintf('L_dpu = %.4f H\n', L_dpu);
L_mdpu = L_md / L_N;
fprintf('L_mdpu = %.4f H\n', L_mdpu);
L_qpu = L_q / L_N;
fprintf('L_qpu = %.4f H\n', L_qpu);
L_mqpu = L_mq / L_N;
fprintf('L_mqpu = %.4f H\n', L_mqpu);
L_ssigmapu = L_ssigma / L_N;
fprintf('L_ssigmapu = %.4f H\n', L_ssigmapu);