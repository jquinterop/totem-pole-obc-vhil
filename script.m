%% =========================================================================
%  RÉPLICA BELKAMEL, KIM Y CHOI (2021)
%  IEEE Trans. Power Electronics, vol. 36, no. 3, pp. 3486-3495.
%
%  Contiene:
%    - PLL SRF-SOGI validado matemáticamente (PI sobre v_d, SOGI fijo en )
%    - Lazo interno de corriente (PI con anti-windup)
%    - Lazo externo de tensión (PI con anti-windup)
%    - Simulación lazo cerrado dinámica (320 V  400 V)
%    - Réplica Fig.12(a) modo G2V y Fig.12(b) modo V2G
%
%  Autor: Juan Pablo Quintero Pachón — Universidad de los Andes
% =========================================================================

clear; clc; close all;

fprintf('=========================================================\n');
fprintf('  RÉPLICA BELKAMEL 2021 — Versión final con PLL SOGI\n');
fprintf('=========================================================\n\n');

%% =========================================================================
%  SECCIÓN 1: PARÁMETROS DEL PAPER (Tabla I + Sección IV)
% =========================================================================

Vg_rms  = 220;
Vg_pk   = Vg_rms * sqrt(2);   % 311.13 V
fg      = 60;
omg     = 2 * pi * fg;         % 376.99 rad/s

P       = 3300;                % W
Vdc_ref = 400;                 % V
Vdc_min = 320;                 % V
Vdc_op  = 350;                 % V — punto de operación

n       = 1.7;
fs      = 100e3;
Ts      = 1 / fs;
Lg      = 350e-6;
Cc      = 2.5e-6;
Ls      = 45e-6;
Cf      = 45e-6;
k_mod   = 0.5;
phi_max =  0.45;
phi_min = -0.45;
Ts_ctrl = 1 / fs;

%% =========================================================================
%  SECCIÓN 2: PUNTO DE OPERACIÓN Y GANANCIAS DE PLANTA (ec.10 paper)
%
%  Ecuación (10) del paper:
%    p(t) = (4²  2 + d²  d + 0.25) · n · Vdc · |vg(t)| / (fs · Ls)
%
%  Con sustitución de signo para _control positivo en G2V:
%    p(t) = (4² + 2 + d²  d + 0.25) · n · Vdc · |vg(t)| / (fs · Ls)
% =========================================================================

theta_ciclo = linspace(0, 2*pi, 100000);
d_ciclo     = k_mod * abs(sin(theta_ciclo));
vg_ciclo    = Vg_pk * abs(sin(theta_ciclo));

P_avg_fn = @(phi) mean( ...
    (4*phi^2 + 2*phi + d_ciclo.^2 - d_ciclo + 0.25) .* ...
     n * Vdc_op .* vg_ciclo / (fs * Ls) );

phi_op = fzero(@(ph) P_avg_fn(ph) - P, 0.1);

fprintf('--- Punto de operación ---\n');
fprintf('  phi_op       = %.4f rad\n', phi_op);
fprintf('  P(phi_op)    = %.0f W (objetivo: %.0f W)\n', P_avg_fn(phi_op), P);

dP_dphi   = (P_avg_fn(phi_op + 1e-5) - P_avg_fn(phi_op - 1e-5)) / (2e-5);
K_i_plant = 2 * dP_dphi / Vg_pk;
Ig_pk_nom = 2 * P / Vg_pk;

fprintf('  dP/dphi      = %.2f W/rad\n', dP_dphi);
fprintf('  K_i_plant    = %.4f A/rad\n', K_i_plant);
fprintf('  Ig_pk_nom    = %.2f A\n\n', Ig_pk_nom);

%% =========================================================================
%  SECCIÓN 3: DISEÑO DEL SRF-PLL SOGI (versión validada)
%
%  ARQUITECTURA:
%    vg  SOGI( fijo)   = y_d,  = y_q
%    v_d = ·cos() + ·sin()  (NO v_q)
%    v_d  PI    +  integrador + wrap(2)  _g
%
%  Ventajas:
%    - Sin offset (_g queda en fase con vg)
%    - Sin normalización (PI ve v_d directo)
%    - SOGI fijo en  (robusto, sin feedback inestable)
%
%  Matemática:
%    Con SOGI convergido:  = Vpk·sin(t),  = -Vpk·cos(t)
%    v_d = Vpk·[sin(t)·cos() - cos(t)·sin()] = Vpk·sin(t-)
%    Para v_d = 0   = t (en fase con vg, sin offset)
%
%  Diseño del PI:
%    H_CL(s) = (Vpk·Kp·s + Vpk·Ki) / (s² + Vpk·Kp·s + Vpk·Ki)
%    Comparando con forma estándar 2° orden: 2n = Vpk·Kp, n² = Vpk·Ki
%    Entonces: Kp = 2n/Vpk, Ki = n²/Vpk
% =========================================================================

% --- SOGI fijo en  ---
k_sogi = 0.75;

% Ecuaciones continuas:
%   G_d(s) = k··s / (s² + k··s + ²)
%   G_q(s) = k·²  / (s² + k··s + ²)
s_tf  = tf('s');
G_d_s = (k_sogi * omg * s_tf) / (s_tf^2 + k_sogi*omg*s_tf + omg^2);
G_q_s = (k_sogi * omg^2)      / (s_tf^2 + k_sogi*omg*s_tf + omg^2);

% Discretización Tustin ( constante  coeficientes fijos)
opts = c2dOptions('Method', 'tustin', 'PrewarpFrequency', omg);
G_d_z = c2d(G_d_s, Ts_ctrl, opts);
G_q_z = c2d(G_q_s, Ts_ctrl, opts);

[SOGI_num_d, SOGI_den_d] = tfdata(G_d_z, 'v');
[SOGI_num_q, SOGI_den_q] = tfdata(G_q_z, 'v');
SOGI_den = SOGI_den_d;  % denominador común

% --- PI del PLL ---
BW_pll      = 30;
zeta_pll    = 0.707;
omega_n_pll = 2 * pi * BW_pll;
Kp_pll      = 2 * zeta_pll * omega_n_pll / Vg_pk;
Ki_pll      = omega_n_pll^2 / Vg_pk;

fprintf('--- PLL SRF-SOGI (fijo en ) ---\n');
fprintf('  k_sogi       = %.2f\n', k_sogi);
fprintf('  SOGI_num_d   = [%.6f, %.6f, %.6f]\n', SOGI_num_d);
fprintf('  SOGI_num_q   = [%.6f, %.6f, %.6f]\n', SOGI_num_q);
fprintf('  SOGI_den     = [%.6f, %.6f, %.6f]\n', SOGI_den);
fprintf('  BW_pll       = %.0f Hz | zeta = %.3f\n', BW_pll, zeta_pll);
fprintf('  Kp_pll       = %.6f\n', Kp_pll);
fprintf('  Ki_pll       = %.6f\n\n', Ki_pll);

%% =========================================================================
%  SECCIÓN 4: DISEÑO DEL LAZO INTERNO DE CORRIENTE
% =========================================================================

BW_i      = 500;
PM_i      = 60;
omega_c_i = 2 * pi * BW_i;
omega_z_i = omega_c_i / tan(PM_i * pi/180);
Kp_i      = 1 / (K_i_plant * sqrt(1 + (omega_z_i/omega_c_i)^2));
Ki_i      = Kp_i * omega_z_i;
Kaw_i     = sqrt(Ki_i / Kp_i);

fprintf('--- PI lazo interno (corriente) ---\n');
fprintf('  BW_i         = %.0f Hz | PM = %d°\n', BW_i, PM_i);
fprintf('  Kp_i         = %.6f\n', Kp_i);
fprintf('  Ki_i         = %.6f\n', Ki_i);
fprintf('  Kaw_i        = %.6f\n\n', Kaw_i);

%% =========================================================================
%  SECCIÓN 5: DISEÑO DEL LAZO EXTERNO DE TENSIÓN
% =========================================================================

I_carga   = P / Vdc_ref;
I_trans   = P_avg_fn(phi_op) / Vdc_op;
I_neta    = I_trans - I_carga;
t_rise    = 0.080;
DeltaVdc  = Vdc_ref - Vdc_min;
C_bat     = max(I_neta * t_rise / DeltaVdc, 1e-3);
C_eff     = C_bat + Cf;

BW_v      = 10;
PM_v      = 60;
omega_c_v = 2 * pi * BW_v;
G_v_wc    = Vg_pk / (2 * Vdc_ref * C_bat * omega_c_v);
omega_z_v = omega_c_v / tan(PM_v * pi/180);
Kp_v      = 1 / (G_v_wc * sqrt(1 + (omega_z_v/omega_c_v)^2));
Ki_v      = Kp_v * omega_z_v;
Kaw_v     = sqrt(Ki_v / Kp_v);
igdref_max = 1.5 * Ig_pk_nom;
igdref_min = 0;

fprintf('--- PI lazo externo (tensión) ---\n');
fprintf('  C_bat        = %.3f mF | C_eff = %.3f mF\n', C_bat*1e3, C_eff*1e3);
fprintf('  BW_v         = %.0f Hz | PM = %d°\n', BW_v, PM_v);
fprintf('  Kp_v         = %.6f\n', Kp_v);
fprintf('  Ki_v         = %.6f\n', Ki_v);
fprintf('  Kaw_v        = %.6f\n\n', Kaw_v);

%% =========================================================================
%  SECCIÓN 6A: VALIDACIÓN AISLADA DEL PLL SRF-SOGI (solo MATLAB)
%
%  Simulamos el PLL aislado con vg ideal para confirmar que:
%    - _g converge en fase con vg (sin offset)
%    - _est converge a omg = 2·60
%
%  Este PLL es el que se llevará tal cual a Simulink.
% =========================================================================

fprintf('--- Verificación aislada del PLL SOGI ---\n');

t_pll_sim = 0.3;
t_pll     = 0:Ts_ctrl:t_pll_sim;
Npll      = length(t_pll);

vg_pll_input = Vg_pk * sin(omg * t_pll);

% Estados del PLL
theta_pll_h = zeros(1, Npll);
omega_est_h = zeros(1, Npll);
theta_pll   = 0;
int_pll     = 0;
yd_p1 = 0; yd_p2 = 0;
yq_p1 = 0; yq_p2 = 0;
vg_p1 = 0; vg_p2 = 0;

for i = 1:Npll
    vg_k = vg_pll_input(i);
    
    % SOGI (dos filtros IIR en paralelo)
    yd_v = -SOGI_den(2)*yd_p1 - SOGI_den(3)*yd_p2 ...
         + SOGI_num_d(1)*vg_k + SOGI_num_d(2)*vg_p1 + SOGI_num_d(3)*vg_p2;
    yq_v = -SOGI_den(2)*yq_p1 - SOGI_den(3)*yq_p2 ...
         + SOGI_num_q(1)*vg_k + SOGI_num_q(2)*vg_p1 + SOGI_num_q(3)*vg_p2;
    
    vg_p2 = vg_p1;   vg_p1 = vg_k;
    yd_p2 = yd_p1;   yd_p1 = yd_v;
    yq_p2 = yq_p1;   yq_p1 = yq_v;
    
    % Park: v_d = ·cos() + ·sin()  SE ALIMENTA AL PI (no v_q)
    vd_pll = yd_v * cos(theta_pll) + yq_v * sin(theta_pll);
    
    % PI del PLL
    int_pll = int_pll + Ki_pll * Ts_ctrl * vd_pll;
    delta_omega = Kp_pll * vd_pll + int_pll;
    omega_est = omg + delta_omega;
    
    % Integrador con wrap 2
    theta_pll = mod(theta_pll + omega_est * Ts_ctrl, 2*pi);
    
    theta_pll_h(i) = theta_pll;
    omega_est_h(i) = omega_est;
end

% Evaluación en régimen permanente
idx_pll_ss = round(0.8*Npll):Npll;
theta_ideal_pll = mod(omg * t_pll, 2*pi);
err_pll = theta_pll_h(idx_pll_ss) - theta_ideal_pll(idx_pll_ss);
err_pll = mod(err_pll + pi, 2*pi) - pi;

fprintf('  Régimen permanente (t > %.2f s):\n', t_pll(idx_pll_ss(1)));
fprintf('    _est = %.4f rad/s (ideal: %.4f)\n', ...
        mean(omega_est_h(idx_pll_ss)), omg);
fprintf('    Error  = %.4f ms (%.4f°)\n\n', ...
        mean(abs(err_pll))/omg*1e3, mean(abs(err_pll))*180/pi);

%% =========================================================================
%  SECCIÓN 6B: SIMULACIÓN DEL LAZO DE CONTROL (MATLAB)
%
%  IMPORTANTE — DECISIÓN DE DISEÑO:
%  En MATLAB usamos  IDEAL para el lazo de control. Esto NO es hacer trampa
%  porque:
%    1. El PLL ya fue validado por separado (Sección 6A) — sabemos que
%       converge correctamente
%    2. En MATLAB queremos validar los PI de tensión y corriente, no el
%       acople transitorio PLL+control
%    3. En Simulink SÍ se usará el PLL real — ahí vemos el sistema completo
%
%  Esta separación permite validar cada subsistema independientemente.
% =========================================================================

fprintf('--- Simulación lazo de control ( ideal, PLL validado aparte) ---\n');

t_sim = 0.100;
dt    = Ts_ctrl;
t_vec = 0:dt:t_sim;
N     = length(t_vec);

% Almacenamiento
Vdc_h    = zeros(1,N);
phi_h    = zeros(1,N);
P_h      = zeros(1,N);
igdref_h = zeros(1,N);
igd_h    = zeros(1,N);

% Estados iniciales
Vdc_s  = Vdc_min;
theta_k = 0;

% Estados del SOGI para extraer i_gd (sobre ig)
yd_i_p1 = 0; yd_i_p2 = 0;
yq_i_p1 = 0; yq_i_p2 = 0;
ig_p1   = 0; ig_p2   = 0;

% Estados de los lazos
int_v    = 0;
igdref_k = 0;
int_i    = 0;
phi_k    = 0;

for i = 1:N
    t_k  = t_vec(i);

    % BLOQUE 1:  ideal (equivalente a PLL ya convergido)
    theta_k = mod(omg * t_k, 2*pi);
    
    % BLOQUE 2: Corriente de red (modelo promediado)
    Ig_pk_act = 2 * max(P_h(max(1,i-1)), 0) / Vg_pk;
    ig_k      = Ig_pk_act * sin(theta_k);
    
    % BLOQUE 3: Extracción de i_gd con SOGI + Park
    yd_i = -SOGI_den(2)*yd_i_p1 - SOGI_den(3)*yd_i_p2 ...
         + SOGI_num_d(1)*ig_k + SOGI_num_d(2)*ig_p1 + SOGI_num_d(3)*ig_p2;
    yq_i = -SOGI_den(2)*yq_i_p1 - SOGI_den(3)*yq_i_p2 ...
         + SOGI_num_q(1)*ig_k + SOGI_num_q(2)*ig_p1 + SOGI_num_q(3)*ig_p2;
    
    ig_p2 = ig_p1;   ig_p1 = ig_k;
    yd_i_p2 = yd_i_p1;  yd_i_p1 = yd_i;
    yq_i_p2 = yq_i_p1;  yq_i_p1 = yq_i;
    
    i_gd_k = yd_i * sin(theta_k) - yq_i * cos(theta_k);
    
    % BLOQUE 4: PI externo — lazo de tensión (con anti-windup)
    err_v    = Vdc_ref - Vdc_s;
    igdref_u = Kp_v * err_v + int_v;
    igdref_k = max(igdref_min, min(igdref_max, igdref_u));
    int_v    = int_v + Ki_v*dt*err_v + Kaw_v*(igdref_k - igdref_u);
    
    % BLOQUE 5: PI interno — lazo de corriente (con anti-windup)
    err_i  = igdref_k - i_gd_k;
    phi_u  = Kp_i * err_i + int_i;
    phi_k  = max(phi_min, min(phi_max, phi_u));
    int_i  = int_i + Ki_i*dt*err_i + Kaw_i*(phi_k - phi_u);
    
    % BLOQUE 6: Ec.(10) del paper
    d_k      = k_mod * abs(sin(theta_k));
    vg_abs_k = Vg_pk * abs(sin(theta_k));
    P_k      = max(0, (4*phi_k^2 + 2*phi_k + d_k^2 - d_k + 0.25) * ...
                   n * Vdc_s * vg_abs_k / (fs * Ls));
    
    % BLOQUE 7: Dinámica de Vdc
    dVdc  = (P_k/max(Vdc_s,1) - I_carga) / C_eff;
    Vdc_s = max(Vdc_min, min(Vdc_ref, Vdc_s + dVdc*dt));
    
    Vdc_h(i)    = Vdc_s;
    phi_h(i)    = phi_k;
    P_h(i)      = P_k;
    igdref_h(i) = igdref_k;
    igd_h(i)    = i_gd_k;
end

idx_ss = round(0.85*N):N;
fprintf('  Régimen permanente (t > %.0f ms):\n', t_vec(idx_ss(1))*1e3);
fprintf('    Vdc  = %.2f V (ref: %.0f V)\n', mean(Vdc_h(idx_ss)), Vdc_ref);
fprintf('    phi  = %.4f rad (phi_op: %.4f rad)\n', mean(phi_h(idx_ss)), phi_op);
fprintf('    P    = %.0f W (nominal: %.0f W)\n\n', mean(P_h(idx_ss)), P);

%% =========================================================================
%  SECCIÓN 7: RÉGIMEN PERMANENTE IDEAL — RÉPLICA FIG.12 DEL PAPER
% =========================================================================

t_ss    = 0 : 1/(2000*fg) : 3/fg;
Ig_pk   = 2 * P / Vg_pk;
ig_G2V  = Ig_pk * sin(omg * t_ss);
vg_ss   = Vg_pk * sin(omg * t_ss);
vCc_ss  = 2 * abs(vg_ss);
idc_G2V = (P/Vdc_ref) * (1 - cos(2*omg*t_ss));
ig_V2G  = -ig_G2V;
idc_V2G = -idc_G2V;

%% =========================================================================
%  SECCIÓN 8: FIGURAS
% =========================================================================

t_ms    = t_vec * 1e3;
t_ss_ms = t_ss  * 1e3;

% ---------- FIGURA 1: Dinámica del lazo de control cascada ----------
figure('Name','Fig.1 — Lazo de control cascada', 'Position',[50 420 1100 720])

subplot(4,1,1)
plot(t_ms, Vdc_h,'b','LineWidth',2); hold on
yline(Vdc_ref,'r--','LineWidth',1.2); yline(Vdc_min,'k:','LineWidth',0.8)
ylabel('V_{dc} (V)'); ylim([310 415])
title('V_{dc}(t) — cargando de 320 V a 400 V'); grid on

subplot(4,1,2)
plot(t_ms, igdref_h,'r','LineWidth',1.5); hold on
plot(t_ms, igd_h,'b','LineWidth',1)
ylabel('i_{gd} (A)')
title('i_{gd}^{*} (referencia) y i_{gd} (medida con SOGI + Park)')
legend('i_{gd}^{*}','i_{gd}','Location','best'); grid on

subplot(4,1,3)
plot(t_ms, phi_h,'m','LineWidth',1.5); hold on
yline(phi_op,'g--',sprintf('\\phi_{op}=%.3f',phi_op),'LineWidth',1)
yline(phi_max,'r:','\phi_{max}','LineWidth',0.8)
ylabel('\phi (rad)'); ylim([-0.02 0.5])
title('\phi(t) — salida del PI interno'); grid on

subplot(4,1,4)
plot(t_ms, P_h/1e3,'k','LineWidth',1.5); hold on
yline(P/1e3,'r--','P_{nom}','LineWidth',1)
xlabel('t (ms)'); ylabel('P (kW)'); ylim([-0.1 4.2])
title('Potencia activa P(t)'); grid on

sgtitle('Fig.1 — Dinámica del lazo de control ( ideal; PLL validado en Fig.4)', ...
        'FontSize',11,'FontWeight','bold')

% ---------- FIGURA 2: Réplica Fig.12(a) — modo G2V ----------
figure('Name','Fig.2 — Régimen permanente G2V', 'Position',[50 50 960 600])

subplot(4,1,1); plot(t_ss_ms,ig_G2V,'r','LineWidth',1.5); hold on
yline(0,'k:'); ylabel('i_g (A)'); ylim([-30 30])
title(sprintf('i_g — I_{pk}=%.1fA — sin pico en cruce por cero',Ig_pk)); grid on

subplot(4,1,2); plot(t_ss_ms,vg_ss,'b','LineWidth',1.5)
ylabel('v_g (V)'); ylim([-400 400]); grid on
title(sprintf('v_g — V_{pk}=%.1fV, f_g=%.0fHz',Vg_pk,fg))

subplot(4,1,3); plot(t_ss_ms,vCc_ss,'m','LineWidth',1.5)
ylabel('v_{Cc} (V)'); ylim([0 750]); grid on
title(sprintf('v_{Cc}=2|v_g| — pico=%.0fV',2*Vg_pk))

subplot(4,1,4); plot(t_ss_ms,idc_G2V,'g','LineWidth',1.5); hold on
yline(P/Vdc_ref,'k--',sprintf('I_{DC}=%.2fA',P/Vdc_ref))
ylabel('i_{dc} (A)'); xlabel('t (ms)'); ylim([-1 20]); grid on
title('i_{dc} — rizo 120 Hz inherente')

sgtitle({'Fig.2 — Réplica Fig.12(a) — Forward mode G2V', ...
         sprintf('V_{dc}=%.0fV | P=%.0fW | \\phi=%.3frad',Vdc_ref,P,phi_op)}, ...
        'FontSize',11,'FontWeight','bold')

% ---------- FIGURA 3: Réplica Fig.12(b) — modo V2G ----------
figure('Name','Fig.3 — Régimen permanente V2G', 'Position',[1020 420 960 600])

subplot(4,1,1); plot(t_ss_ms,ig_V2G,'r','LineWidth',1.5); hold on
yline(0,'k:'); ylabel('i_g (A)'); ylim([-30 30])
title('i_g — corriente invertida (batería inyecta a red)'); grid on

subplot(4,1,2); plot(t_ss_ms,vg_ss,'b','LineWidth',1.5)
ylabel('v_g (V)'); ylim([-400 400]); grid on; title('v_g')

subplot(4,1,3); plot(t_ss_ms,vCc_ss,'m','LineWidth',1.5)
ylabel('v_{Cc} (V)'); ylim([0 750]); grid on; title('v_{Cc}=2|v_g|')

subplot(4,1,4); plot(t_ss_ms,idc_V2G,'g','LineWidth',1.5); hold on
yline(-P/Vdc_ref,'k--',sprintf('I_{DC}=%.2fA',P/Vdc_ref))
ylabel('i_{dc} (A)'); xlabel('t (ms)'); ylim([-20 1]); grid on
title('i_{dc} — flujo inverso')

sgtitle({'Fig.3 — Réplica Fig.12(b) — Reverse mode V2G', ...
         sprintf('V_{dc}=%.0fV | P=%.0fW | \\phi=%.3frad',Vdc_ref,P,-phi_op)}, ...
        'FontSize',11,'FontWeight','bold')

% ---------- FIGURA 4: Verificación del PLL SOGI (aislado, Sección 6A) ----------
theta_ideal_pll_full = mod(omg * t_pll, 2*pi);
t_pll_ms = t_pll * 1e3;

figure('Name','Fig.4 — Verificación del PLL SOGI', 'Position',[550 100 1100 700])

subplot(3,1,1)
plot(t_pll_ms, vg_pll_input,'b','LineWidth',1.2); hold on
ylabel('v_g (V)'); ylim([-400 400])
title('Tensión de red v_g (entrada al PLL)'); grid on
xlim([0 100])

subplot(3,1,2)
plot(t_pll_ms, theta_pll_h,'b','LineWidth',1.5); hold on
plot(t_pll_ms, theta_ideal_pll_full,'r--','LineWidth',1)
ylabel('\theta (rad)'); ylim([0 2*pi])
title('_g PLL (azul) vs  ideal (rojo) — deben superponerse tras transitorio')
legend('\theta_g PLL','\theta ideal','Location','best'); grid on

subplot(3,1,3)
err_all = theta_pll_h - theta_ideal_pll_full;
err_all = mod(err_all + pi, 2*pi) - pi;
plot(t_pll_ms, err_all*180/pi,'k','LineWidth',1.2); hold on
yline(0,'r:'); xlabel('t (ms)'); ylabel('Error  (°)')
title(sprintf('Error de seguimiento del PLL — error medio en régimen: %.3f°', ...
              mean(abs(err_all(idx_pll_ss)))*180/pi))
grid on

sgtitle('Fig.4 — Validación aislada del PLL SRF-SOGI (se llevará a Simulink)', ...
        'FontSize',11,'FontWeight','bold')

%% =========================================================================
%  SECCIÓN 9: EXPORTAR PARÁMETROS PARA SIMULINK
% =========================================================================

bp = struct();

% --- Red ---
bp.Vg_rms = Vg_rms;  bp.Vg_pk = Vg_pk;
bp.fg = fg;  bp.omg = omg;

% --- Potencia ---
bp.P_rated = P;
bp.Vdc_ref = Vdc_ref;  bp.Vdc_min = Vdc_min;  bp.Vdc_op = Vdc_op;

% --- Convertidor ---
bp.n  = n;   bp.fs = fs;   bp.Ts = Ts;
bp.Lg = Lg;  bp.Cc = Cc;   bp.Ls = Ls;   bp.Cf = Cf;
bp.k_mod   = k_mod;
bp.phi_max = phi_max;  bp.phi_min = phi_min;  bp.phi_op = phi_op;

% --- Transformador ---
bp.L1_trafo  = n^2 * Ls;
bp.M_trafo   = 0.99 * n * Ls;
bp.L2_trafo  = Ls;
bp.Rc_trafo1 = 0.1;  bp.Rc_trafo2 = 0.06;

% --- Semiconductores ---
bp.Ron_mos = 0.045;  bp.Vf_body = 1.2;
bp.Rg = 0.03;

% --- PLL SRF-SOGI ---
bp.k_sogi     = k_sogi;
bp.SOGI_num_d = SOGI_num_d;
bp.SOGI_num_q = SOGI_num_q;
bp.SOGI_den   = SOGI_den;
bp.BW_pll     = BW_pll;
bp.Kp_pll     = Kp_pll;
bp.Ki_pll     = Ki_pll;

% --- Lazos de control ---
bp.Kp_v = Kp_v;  bp.Ki_v = Ki_v;  bp.Kaw_v = Kaw_v;
bp.igdref_max = igdref_max;  bp.igdref_min = igdref_min;
bp.Kp_i = Kp_i;  bp.Ki_i = Ki_i;  bp.Kaw_i = Kaw_i;
bp.Ts_ctrl = Ts_ctrl;

% --- Carga ---
bp.I_load = P / Vdc_ref;
bp.C_bat  = C_bat;  bp.C_eff = C_eff;

assignin('base', 'belkamel_params', bp);

sep = repmat('=', 1, 60);
fprintf('%s\n  PARÁMETROS EXPORTADOS A belkamel_params\n%s\n\n', sep, sep);

fprintf('--- Red ---\n');
fprintf('  Vg_rms    = %.1f Vrms  | Vg_pk = %.2f V\n', bp.Vg_rms, bp.Vg_pk);
fprintf('  fg        = %.0f Hz      | omg   = %.2f rad/s\n', bp.fg, bp.omg);

fprintf('\n--- Potencia ---\n');
fprintf('  P_rated   = %.0f W\n', bp.P_rated);
fprintf('  Vdc_ref   = %.0f V      | Vdc_min = %.0f V\n', bp.Vdc_ref, bp.Vdc_min);

fprintf('\n--- Convertidor ---\n');
fprintf('  n = %.1f | fs = %.0f kHz | Ts = %.1f us\n', bp.n, bp.fs/1e3, bp.Ts*1e6);
fprintf('  Lg = %.0f uH | Cc = %.1f uF | Ls = %.0f uH | Cf = %.0f uF\n', ...
        bp.Lg*1e6, bp.Cc*1e6, bp.Ls*1e6, bp.Cf*1e6);
fprintf('  phi_op    = %.4f rad\n', bp.phi_op);

fprintf('\n--- PLL SRF-SOGI ---\n');
fprintf('  k_sogi    = %.2f\n', bp.k_sogi);
fprintf('  SOGI_num_d = [%.6f, %.6f, %.6f]\n', bp.SOGI_num_d);
fprintf('  SOGI_num_q = [%.6f, %.6f, %.6f]\n', bp.SOGI_num_q);
fprintf('  SOGI_den   = [%.6f, %.6f, %.6f]\n', bp.SOGI_den);
fprintf('  Kp_pll    = %.6f | Ki_pll = %.6f\n', bp.Kp_pll, bp.Ki_pll);

fprintf('\n--- Lazos de control ---\n');
fprintf('  [Tensión] Kp_v = %.6f | Ki_v = %.6f | Kaw_v = %.6f\n', ...
        bp.Kp_v, bp.Ki_v, bp.Kaw_v);
fprintf('  [Corriente] Kp_i = %.6f | Ki_i = %.6f | Kaw_i = %.6f\n', ...
        bp.Kp_i, bp.Ki_i, bp.Kaw_i);
fprintf('  Ts_ctrl   = %.1f us\n', bp.Ts_ctrl*1e6);

fprintf('\n--- Carga ---\n');
fprintf('  I_load    = %.3f A   | C_bat = %.3f mF\n', bp.I_load, bp.C_bat*1e3);

fprintf('\n%s\n\n', sep);
fprintf('  INSTRUCCIONES PARA SIMULINK:\n');
fprintf('  ============================\n\n');
fprintf('  [Subsistema SRF_PLL_SOGI]\n');
fprintf('    Dos bloques "Discrete Transfer Fcn":\n');
fprintf('      SOGI y_d:  num = belkamel_params.SOGI_num_d\n');
fprintf('                 den = belkamel_params.SOGI_den\n');
fprintf('      SOGI y_q:  num = belkamel_params.SOGI_num_q\n');
fprintf('                 den = belkamel_params.SOGI_den\n');
fprintf('      Sample time = belkamel_params.Ts_ctrl en ambos\n\n');
fprintf('    Park dq: v_d = ·cos() + ·sin()   ALIMENTAR PI AL V_D\n');
fprintf('    PI: Kp = belkamel_params.Kp_pll, Ki = belkamel_params.Ki_pll\n');
fprintf('    Suma + omg = belkamel_params.omg\n');
fprintf('    Integrador Forward Euler + mod(2*pi)\n');
fprintf('    Realimentación sin() y cos() al Park\n\n');
fprintf('  [Subsistema Extraccion_igd]\n');
fprintf('    Misma estructura de SOGI, pero sobre ig\n');
fprintf('    Park: i_gd = ·cos(_g) + ·sin(_g), _g viene del PLL\n\n');
fprintf('%s\n', sep);