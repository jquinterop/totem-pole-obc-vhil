%% =========================================================================
%  CONVERTIDOR OBC BIDIRECCIONAL ISOLADO — TESIS IELE3002
%  Caso: 120 Vac / 60 Hz, 1 kW, Vdc = 300-450 V
%
%  Topología base: Belkamel, Kim y Choi (2021)
%  IEEE Trans. Power Electronics, vol. 36, no. 3, pp. 3486-3495.
%
%  Autor: Juan Pablo Quintero Pachón — Universidad de los Andes
% =========================================================================

clear; clc; close all;

fprintf('=========================================================\n');
fprintf('  CONVERTIDOR OBC IELE3002 — 120V/60Hz, 1kW\n');
fprintf('=========================================================\n\n');

sat = @(x,xmin,xmax) min(max(x,xmin),xmax);

%% =========================================================================
%  SECCIÓN 1: ESPECIFICACIONES (IFEC 2026)
% =========================================================================

% --- Lado AC (red) ---
Vg_rms  = 120;                       % IFEC: 120 Vac
Vg_pk   = Vg_rms * sqrt(2);          % = 169.7 V
fg      = 60;                        % IFEC: 60 Hz
omg     = 2*pi*fg;

% --- Lado DC (batería) ---
P        = 1000;                     % IFEC: 1 kW máximo
Vdc_ref  = 400;                      % Centro del rango operativo
Vdc_min  = 300;                      % IFEC: límite inferior
Vdc_max  = 450;                      % IFEC: límite superior
Vbat_nom = 380;                      % Batería simulada (modificable)

% --- Switching ---
fs      = 100e3;                     % 100 kHz (SiC práctico)
Ts      = 1/fs;
Ts_ctrl = Ts;                        % Control sincrónico con switching

% --- Modulación SPWM (Belkamel) ---
k_mod   = 0.5;                       % Constante del paper, no se modifica

fprintf('--- Especificaciones IFEC ---\n');
fprintf('  Vg_rms = %.0f V | fg = %.0f Hz | Vg_pk = %.1f V\n', Vg_rms, fg, Vg_pk);
fprintf('  P_rated = %.0f W | Vdc = %.0f-%.0f V (ref %.0f V)\n', ...
        P, Vdc_min, Vdc_max, Vdc_ref);
fprintf('  fs = %.0f kHz\n\n', fs/1e3);

%% =========================================================================
%  SECCIÓN 2: DIMENSIONAMIENTO ANALÍTICO DE COMPONENTES
% =========================================================================

fprintf('=========================================================\n');
fprintf('  DIMENSIONAMIENTO COMPONENTE POR COMPONENTE\n');
fprintf('=========================================================\n\n');

% -------------------------------------------------------------------------
% 2.1 Relación de vueltas del transformador HF (n = primario:secundario)
%
% Criterio: n·Vdc_max <= vCc_pk = 2·Vg_pk
% Para evitar saturar el modulador cuando vCc esta en su minimo.
% En 120V, vCc_pk = 340V, Vdc_max = 450V => n_max ~ 0.756
% -------------------------------------------------------------------------

vCc_pk = 2*Vg_pk;                              % = 339.4 V
n_max  = vCc_pk / Vdc_max;
n      = 0.7;                                  % Margen ~8% sobre n_max
margin_n = (n_max - n)/n_max * 100;

fprintf('--- 2.1 Relación de vueltas n ---\n');
fprintf('  vCc_pk         = 2·Vg_pk = %.1f V\n', vCc_pk);
fprintf('  n_max teorico  = vCc_pk/Vdc_max = %.3f\n', n_max);
fprintf('  n adoptado     = %.2f (margen %.1f%%)\n', n, margin_n);
fprintf('  ATENCION: n < 1 => trafo "step-up" (inverso al paper donde n=1.7)\n\n');

% -------------------------------------------------------------------------
% 2.2 Inductores totem-pole interleaved Lg1, Lg2
%
% Criterio: ripple de corriente <= 12% de I_g_pk para preservar ZVS
% en todo el rango. Formula del ripple en boost interleaved con
% 50% duty fijo:
%   dI_Lg = (Vg_pk · 0.5) / (Lg · fs · 2)
% (factor 1/2 por la simetria del switching en 50% duty)
% -------------------------------------------------------------------------

Ig_pk     = 2*P/Vg_pk;                         % = 11.78 A
dI_target = 0.12 * Ig_pk;                      % objetivo: <12% ripple
Lg_min    = (Vg_pk * 0.5) / (dI_target * fs * 2);
Lg        = 300e-6;                            % redondeado al estandar
dI_Lg     = (Vg_pk * 0.5) / (Lg * fs * 2);
dI_Lg_pct = dI_Lg/Ig_pk * 100;

fprintf('--- 2.2 Inductores Lg1, Lg2 ---\n');
fprintf('  Ig_pk          = %.2f A\n', Ig_pk);
fprintf('  Lg_min teorico = %.0f uH (ripple 12%%)\n', Lg_min*1e6);
fprintf('  Lg adoptado    = %.0f uH\n', Lg*1e6);
fprintf('  Ripple real    = %.2f A (%.1f%% de Ig_pk)\n\n', dI_Lg, dI_Lg_pct);

% -------------------------------------------------------------------------
% 2.3 Capacitor de clamp Cc
%
% Criterio: la constante de tiempo Cc·R_carga_eq << 1/(2·fg)
% Para que vCc siga el doble del rectificado |vg| sin atraso significativo.
% Escalado desde el paper: Cc ∝ P / (Vg_pk² · fg)
%
% Comparacion con paper (220V, 3.3kW, fg=60Hz):
%   Cc_paper = 2.5 uF
%   Factor: (P_nuevo/P_paper) · (Vg_pk_paper/Vg_pk_nuevo)^2
% -------------------------------------------------------------------------

Cc_paper      = 2.5e-6;
P_paper       = 3300;
Vgpk_paper    = 220*sqrt(2);
scale_Cc      = (P/P_paper) * (Vgpk_paper/Vg_pk)^2;
Cc_calc       = Cc_paper * scale_Cc;
Cc            = 2.5e-6;                        % coincide casi exacto

fprintf('--- 2.3 Capacitor de clamp Cc ---\n');
fprintf('  Factor escalado = %.3f (P y V cancelan casi exacto)\n', scale_Cc);
fprintf('  Cc calculado    = %.2f uF\n', Cc_calc*1e6);
fprintf('  Cc adoptado     = %.1f uF (valor estandar)\n\n', Cc*1e6);

% -------------------------------------------------------------------------
% 2.4 Inductor serie Ls (transferencia de potencia)
%
% Criterio: ubicar phi_op nominal en el centro del rango util (~0.4 rad)
% para tener margen tanto para regulacion como para perdidas no modeladas.
%
% Formula simplificada (Eq. 10 paper, evaluada en el pico de vg):
%   P_max ≈ (Vab_pk · Vcd_pk · phi) / (2·pi · fs · Ls)
%
% Donde Vab_pk = 2·Vg_pk, Vcd_pk = n·Vdc
% -------------------------------------------------------------------------

phi_op_target = 0.4;                           % rad, centro del rango util
Vab_pk        = 2*Vg_pk;                       % = 340 V
Vcd_pk        = n*Vdc_ref;                     % = 280 V
Ls_calc       = (Vab_pk * Vcd_pk * phi_op_target) / (2*pi*fs*P);
Ls            = 60e-6;                         % redondeado

fprintf('--- 2.4 Inductor serie Ls ---\n');
fprintf('  Vab_pk           = %.1f V\n', Vab_pk);
fprintf('  Vcd_pk = n·Vdc   = %.1f V\n', Vcd_pk);
fprintf('  phi_op objetivo  = %.2f rad\n', phi_op_target);
fprintf('  Ls calculado     = %.1f uH\n', Ls_calc*1e6);
fprintf('  Ls adoptado      = %.0f uH\n\n', Ls*1e6);

% -------------------------------------------------------------------------
% 2.5 Capacitor del bus DC Cf
%
% Criterio: rizado de Vdc a 2·fg (120 Hz) acotado a ≤5%.
% IFEC no especifica rizado, pero 5% es razonable para battery charging.
%
%   dVdc = P / (2·omega·Vdc·Cf)
%   Cf >= P / (2·omega·Vdc²·dVdc_pct)
% -------------------------------------------------------------------------

dVdc_pct_target = 0.05;
Cf_min          = P / (2*omg * Vdc_ref^2 * dVdc_pct_target);
Cf              = 220e-6;                      % valor estandar
dVdc_real       = P / (2*omg * Vdc_ref * Cf);
dVdc_pct_real   = dVdc_real / Vdc_ref * 100;

fprintf('--- 2.5 Capacitor del bus DC Cf ---\n');
fprintf('  Objetivo rizado = %.1f%% de Vdc\n', dVdc_pct_target*100);
fprintf('  Cf_min calc.    = %.0f uF\n', Cf_min*1e6);
fprintf('  Cf adoptado     = %.0f uF (valor estandar)\n', Cf*1e6);
fprintf('  Rizado real     = %.2f V (%.2f%%)\n\n', dVdc_real, dVdc_pct_real);

% -------------------------------------------------------------------------
% 2.6 C_bus de simulacion (estabilizacion numerica Simulink)
%
% Tipicamente 10x el Cf real para estabilizar transitorios durante el
% desarrollo. Cuando el modelo este maduro, reducir hacia Cf real.
% -------------------------------------------------------------------------

C_bus     = 2200e-6;                           % 10x el Cf real
L_par     = 5e-6;                              % parasitica
R_ESR     = 5e-3;                              % ESR cap
f_res     = 1/(2*pi*sqrt(L_par*C_bus));

fprintf('--- 2.6 Filtro de bus para simulacion ---\n');
fprintf('  C_bus           = %.0f uF (10x Cf, estabilizacion numerica)\n', C_bus*1e6);
fprintf('  L_par           = %.1f uH\n', L_par*1e6);
fprintf('  R_ESR           = %.1f mOhm\n', R_ESR*1e3);
fprintf('  f_res LC        = %.1f Hz\n\n', f_res);

C_eff = C_bus;

%% =========================================================================
%  SECCIÓN 3: MODELO PROMEDIADO DE POTENCIA
%
%  Ecuacion del paper (10), con convencion de signos adaptada a Simulink:
%    p(t) = (4·phi² + 2·phi + d² - d + 0.25)·n·Vdc·|vg|/(fs·Ls)
%
%  Calibracion: cuando se ejecute la simulacion en Simulink, ajustar
%  Ls_loss_factor hasta que phi_op del modelo coincida con phi_op de
%  Simulink (depende de perdidas reales, parasitas, etc.).
% =========================================================================

% Factor de calibracion empirica (ajustar tras corrida en Simulink)
Ls_loss_factor = 1;                            % arrancar en 1, ajustar luego
Ls_eff         = Ls * Ls_loss_factor;

theta_ciclo = linspace(0, 2*pi, 200000);
d_ciclo     = k_mod * abs(sin(theta_ciclo));
vg_ciclo    = Vg_pk * abs(sin(theta_ciclo));

S_vg   = mean(vg_ciclo);
S_base = mean((d_ciclo.^2 - d_ciclo + 0.25).*vg_ciclo);

P_avg_model = @(phi,Vdc) ...
    ((4.*phi.^2 + 2.*phi).*S_vg + S_base) .* n .* Vdc ./ (fs.*Ls_eff);

% Punto de operacion calculado
phi_op = fzero(@(ph) P_avg_model(ph,Vdc_ref) - P, [0 2]);

% Linealizacion local
dP_dphi   = (P_avg_model(phi_op + 1e-5,Vdc_ref) - ...
             P_avg_model(phi_op - 1e-5,Vdc_ref)) / (2e-5);
K_i_plant = 2*dP_dphi/Vg_pk;

% Saturacion de phi (con margen para transitorios y perdidas no modeladas)
phi_max = 1.5;                                 % G2V + V2G
phi_min = -1.5;

fprintf('=========================================================\n');
fprintf('  MODELO PROMEDIADO Y PUNTO DE OPERACION\n');
fprintf('=========================================================\n\n');

fprintf('--- Modelo de potencia ---\n');
fprintf('  Ls_loss_factor = %.1f (ajustar tras simulacion en Simulink)\n', Ls_loss_factor);
fprintf('  Ls_eff         = %.1f uH\n', Ls_eff*1e6);
fprintf('  phi_op @ 400V  = %.4f rad\n', phi_op);
fprintf('  P(phi_op)      = %.1f W\n', P_avg_model(phi_op,Vdc_ref));
fprintf('  dP/dphi        = %.2f W/rad\n', dP_dphi);
fprintf('  K_i_plant      = %.4f A/rad\n', K_i_plant);
fprintf('  phi_max/min    = ±%.2f rad\n\n', phi_max);

%% =========================================================================
%  SECCIÓN 4: PLL SRF-SOGI
% =========================================================================

k_sogi     = 0.75;
omega_SOGI = k_sogi * omg / 2;
BW_SOGI    = omega_SOGI/(2*pi);

s_tf  = tf('s');
G_d_s = (k_sogi*omg*s_tf)/(s_tf^2 + k_sogi*omg*s_tf + omg^2);
G_q_s = (k_sogi*omg^2)/(s_tf^2 + k_sogi*omg*s_tf + omg^2);

opts  = c2dOptions('Method','tustin','PrewarpFrequency',omg);
G_d_z = c2d(G_d_s, Ts_ctrl, opts);
G_q_z = c2d(G_q_s, Ts_ctrl, opts);

[SOGI_num_d, SOGI_den_d] = tfdata(G_d_z,'v');
[SOGI_num_q, SOGI_den_q] = tfdata(G_q_z,'v');
SOGI_den = SOGI_den_d;

BW_pll      = 30;
zeta_pll    = 0.707;
omega_n_pll = 2*pi*BW_pll;

Kp_pll = 2*zeta_pll*omega_n_pll/Vg_pk;
Ki_pll = omega_n_pll^2/Vg_pk;

fprintf('--- PLL SRF-SOGI ---\n');
fprintf('  k_sogi  = %.2f | BW_SOGI = %.1f Hz\n', k_sogi, BW_SOGI);
fprintf('  BW_pll  = %.1f Hz\n', BW_pll);
fprintf('  Kp_pll  = %.6f | Ki_pll = %.6f\n\n', Kp_pll, Ki_pll);

%% =========================================================================
%  SECCIÓN 5: DISEÑO DEL LAZO INTERNO DE CORRIENTE
%
%  BW_i = 200 Hz (limitado por BW_SOGI ~70 Hz, se acepta consumo de PM)
%  PM_i = 60 grados
% =========================================================================

BW_i      = 200;
PM_i      = 60;
omega_c_i = 2*pi*BW_i;
omega_z_i = omega_c_i/tan(PM_i*pi/180);

Kp_i  = 1/(K_i_plant*sqrt(1 + (omega_z_i/omega_c_i)^2));
Ki_i  = Kp_i*omega_z_i;
Kaw_i = 1/Kp_i;

fprintf('--- PI lazo interno (corriente) ---\n');
fprintf('  BW_i = %.0f Hz | PM_i = %.0f deg\n', BW_i, PM_i);
fprintf('  Kp_i = %.6f | Ki_i = %.6f | Kaw_i = %.6f\n\n', Kp_i, Ki_i, Kaw_i);

%% =========================================================================
%  SECCIÓN 6: DISEÑO DEL LAZO EXTERNO DE TENSIÓN
%
%  Planta: integrador puro G_v(s) = Vg_pk/(2·Vdc·C_bus·s)
%  BW_v = 5 Hz (bien por debajo de 120 Hz para rechazar rizado de red)
% =========================================================================

BW_v      = 5;
PM_v      = 60;
omega_c_v = 2*pi*BW_v;
omega_z_v = omega_c_v/tan(PM_v*pi/180);

G_v_wc = Vg_pk/(2*Vdc_ref*C_eff*omega_c_v);

Kp_v  = 1/(G_v_wc*sqrt(1 + (omega_z_v/omega_c_v)^2));
Ki_v  = Kp_v*omega_z_v;
Kaw_v = 1/Kp_v;

igdref_max =  1.5*Ig_pk;
igdref_min = -1.5*Ig_pk;

fprintf('--- PI lazo externo (tension) ---\n');
fprintf('  BW_v = %.0f Hz | PM_v = %.0f deg\n', BW_v, PM_v);
fprintf('  Kp_v = %.6f | Ki_v = %.6f | Kaw_v = %.6f\n', Kp_v, Ki_v, Kaw_v);
fprintf('  igdref_max/min = ±%.2f A (simetrica para V2G)\n\n', igdref_max);

%% =========================================================================
%  SECCIÓN 7: SIMULACION DEL LAZO PROMEDIADO (validacion del control)
% =========================================================================

fprintf('=========================================================\n');
fprintf('  SIMULACION DEL LAZO PROMEDIADO\n');
fprintf('=========================================================\n\n');

t_sim = 0.5;
dt    = Ts_ctrl;
t_vec = 0:dt:t_sim;
N     = length(t_vec);

Vdc_h     = zeros(1,N);
phi_h     = zeros(1,N);
P_h       = zeros(1,N);
igdref_h  = zeros(1,N);
igd_h     = zeros(1,N);

Vdc_s = Vdc_min;
int_v = 0; int_i = 0;
igd_s = 0; phi_k = 0; igd_meas = 0;

tau_i = 1/(2*pi*BW_i);
tau_sogi = 1/omega_SOGI;
I_carga  = P/Vdc_ref;

for i = 1:N
    err_v    = Vdc_ref - Vdc_s;
    igdref_u = Kp_v*err_v + int_v;
    igdref_k = sat(igdref_u, igdref_min, igdref_max);
    int_v    = int_v + dt*(Ki_v*err_v + Kaw_v*(igdref_k - igdref_u));

    err_i = igdref_k - igd_meas;
    phi_u = Kp_i*err_i + int_i;
    phi_k = sat(phi_u, phi_min, phi_max);
    int_i = int_i + dt*(Ki_i*err_i + Kaw_i*(phi_k - phi_u));

    P_k    = P_avg_model(phi_k, Vdc_s);
    igd_eq = 2*P_k/Vg_pk;
    igd_s  = igd_s + dt*(igd_eq - igd_s)/tau_i;
    igd_meas = igd_meas + dt*(igd_s - igd_meas)/tau_sogi;

    dVdc  = (P_k/max(Vdc_s,1) - I_carga)/C_eff;
    Vdc_s = sat(Vdc_s + dVdc*dt, 250, 500);

    Vdc_h(i)    = Vdc_s;
    phi_h(i)    = phi_k;
    P_h(i)      = P_k;
    igdref_h(i) = igdref_k;
    igd_h(i)    = igd_meas;
end

idx_ss = round(0.85*N):N;
fprintf('Regimen permanente t > %.0f ms:\n', t_vec(idx_ss(1))*1e3);
fprintf('  Vdc medio = %.2f V (ref = %.0f V, error = %.2f V)\n', ...
        mean(Vdc_h(idx_ss)), Vdc_ref, Vdc_ref - mean(Vdc_h(idx_ss)));
fprintf('  phi medio = %.4f rad\n', mean(phi_h(idx_ss)));
fprintf('  P medio   = %.1f W (nominal %.0f W)\n\n', mean(P_h(idx_ss)), P);

%% =========================================================================
%  SECCIÓN 8: FIGURAS
% =========================================================================

t_ms = t_vec*1e3;

figure('Name','Validacion del control','Position',[100 100 1100 700])

subplot(4,1,1)
plot(t_ms,Vdc_h,'b','LineWidth',2); hold on
yline(Vdc_ref,'r--','V_{dc,ref}','LineWidth',1.2)
yline(Vdc_min,'k:','V_{dc,min}','LineWidth',0.8)
yline(Vdc_max,'k:','V_{dc,max}','LineWidth',0.8)
ylabel('V_{dc} (V)')
ylim([280 470])
title('V_{dc}(t)')
grid on

subplot(4,1,2)
plot(t_ms,igdref_h,'r','LineWidth',1.5); hold on
plot(t_ms,igd_h,'b','LineWidth',1.2)
yline(Ig_pk,'g--',sprintf('I_{g,pk}=%.1f A',Ig_pk),'LineWidth',1)
ylabel('i_{gd} (A)')
title('i_{gd}^{*} (ref) e i_{gd} (medida)')
legend('i_{gd}^{*}','i_{gd}','Location','best')
grid on

subplot(4,1,3)
plot(t_ms,phi_h,'m','LineWidth',1.5); hold on
yline(phi_op,'g--',sprintf('\\phi_{op}=%.3f',phi_op),'LineWidth',1)
yline(phi_max,'r:','\phi_{max}','LineWidth',0.8)
ylabel('\phi (rad)')
ylim([-0.1 1.6])
title('\phi(t)')
grid on

subplot(4,1,4)
plot(t_ms,P_h,'k','LineWidth',1.5); hold on
yline(P,'r--','P_{nom}','LineWidth',1)
xlabel('t (ms)')
ylabel('P (W)')
ylim([0 1300])
title('Potencia transferida')
grid on

sgtitle('Validacion lazo de control IFEC 120V/1kW','FontSize',11,'FontWeight','bold')

%% =========================================================================
%  SECCIÓN 9: EXPORTAR PARAMETROS A WORKSPACE
% =========================================================================

bp = struct();

% Red
bp.Vg_rms = Vg_rms;
bp.Vg_pk  = Vg_pk;
bp.fg     = fg;
bp.omg    = omg;

% Potencia
bp.P_rated = P;
bp.Vdc_ref = Vdc_ref;
bp.Vdc_min = Vdc_min;
bp.Vdc_max = Vdc_max;

% Bateria
bp.Vbat_nom  = Vbat_nom;
bp.Bat_R_int = 0.05;
bp.Bat_SOC0  = 100;
bp.Bat_Q     = 50;
bp.Bat_tau   = 1e-3;

% Convertidor (NUEVOS VALORES para IFEC 120V/1kW)
bp.n  = n;
bp.fs = fs;
bp.Ts = Ts;

bp.Lg = Lg;          % 300 uH
bp.Cc = Cc;          % 2.5 uF
bp.Ls = Ls;          % 60 uH
bp.Cf = Cf;          % 220 uF
bp.Ls_eff = Ls_eff;
bp.Ls_loss_factor = Ls_loss_factor;

bp.k_mod = k_mod;
bp.phi_max = phi_max;
bp.phi_min = phi_min;
bp.phi_op  = phi_op;
bp.K_i_plant = K_i_plant;
bp.Ig_pk_nom = Ig_pk;

% Filtro bus DC
bp.C_bus  = C_bus;
bp.L_par  = L_par;
bp.R_ESR  = R_ESR;
bp.f_res  = f_res;
bp.DeltaVdc     = dVdc_real;
bp.DeltaVdc_pct = dVdc_pct_real;

% Transformador (calculado para n=0.7)
bp.L1_trafo  = n^2*Ls;
bp.M_trafo   = 0.99*n*Ls;
bp.L2_trafo  = Ls;
bp.Rc_trafo1 = 0.1;
bp.Rc_trafo2 = 0.06;

% Semiconductores
bp.Ron_mos = 0.045;
bp.Vf_body = 1.2;
bp.Rg      = 0.03;

% PLL
bp.k_sogi     = k_sogi;
bp.SOGI_num_d = SOGI_num_d;
bp.SOGI_num_q = SOGI_num_q;
bp.SOGI_den   = SOGI_den;
bp.BW_SOGI    = BW_SOGI;
bp.BW_pll     = BW_pll;
bp.Kp_pll     = Kp_pll;
bp.Ki_pll     = Ki_pll;

% Lazos
bp.BW_v = BW_v;
bp.BW_i = BW_i;
bp.Kp_v  = Kp_v;  bp.Ki_v  = Ki_v;  bp.Kaw_v = Kaw_v;
bp.Kp_i  = Kp_i;  bp.Ki_i  = Ki_i;  bp.Kaw_i = Kaw_i;
bp.igdref_max = igdref_max;
bp.igdref_min = igdref_min;
bp.Ts_ctrl = Ts_ctrl;

% Carga
bp.I_load = I_carga;
bp.C_eff  = C_eff;

assignin('base','belkamel_params',bp);

fprintf('=========================================================\n');
fprintf('  PARAMETROS EXPORTADOS A belkamel_params (IFEC 120V/1kW)\n');
fprintf('=========================================================\n\n');

fprintf('Resumen de cambios respecto al paper original (220V/3.3kW):\n');
fprintf('  Vg_rms:   220 -> 120 V\n');
fprintf('  P:       3300 -> 1000 W\n');
fprintf('  n:        1.7 -> %.2f  (TRAFO STEP-UP, invertido)\n', n);
fprintf('  Lg:       350 -> %.0f uH\n', Lg*1e6);
fprintf('  Cc:       2.5 -> %.1f uF  (igual por coincidencia de escalado)\n', Cc*1e6);
fprintf('  Ls:        45 -> %.0f uH\n', Ls*1e6);
fprintf('  Cf:        45 -> %.0f uF  (mayor, IFEC no es capacitor-less)\n', Cf*1e6);
fprintf('  Vdc_ref:  350 -> 400 V\n');
fprintf('\n');

fprintf('SIGUIENTE: simular en Simulink con estos valores. Si phi_op\n');
fprintf('real difiere del calculado (%.3f rad), ajustar Ls_loss_factor\n', phi_op);
fprintf('en Seccion 3 hasta que coincidan.\n\n');