%% Proyecto Global de Asignatura
%% Rubén A. y David A. - Otoño 2024

clear;
decimales = 4;
muestras_para_media = 50;
aproximacion = -3;
%% Fase 1 - Caracterización de la planta
% 
% 
% *Ejercicio 1*
% 
% Capture y caracterice la salida de la planta en función de la tensión de entrada 
% con un  procedimiento análogo al empleado en la Práctica 2 y obtenga su función 
% de transferencia  (FdeT) en las unidades que considere más adecuadas para trabajar 
% en el desarrollo de todo  el diseño. 

ficheros = dir("Fase1\sin_perturbacion\");

funciones_de_transferencia = {};
Kms  = {};
taus = {};
tiempos_de_interes = {};
salidas_de_interes = {};
for i = 1:length(ficheros)
    fichero = ficheros(i).name;
    
    % Saltar '.' y '..'
    if strcmp(fichero, '.') || strcmp(fichero, '..')
        continue;
    end
    str_amplitud_escalon = strsplit(fichero, '_'); str_amplitud_escalon = str_amplitud_escalon{end};
    amplitud_escalon = str2double(regexprep(str_amplitud_escalon, '[^-\d]', ''));

    escalon_sin_perturbacion = load("Fase1\sin_perturbacion\" + fichero);
    Ts = 35e-3;
    
    tiempo     = escalon_sin_perturbacion(:, 1);
    referencia = escalon_sin_perturbacion(:, 2);
    salida     = escalon_sin_perturbacion(:, 3);
    
    x_ini = find(abs(referencia) > 0, 1, 'first') + 1; % Se suma 1 ya que se tienen dos retardos y se va a trabajar con una aproximación de primer orden
    x_fin = length(referencia);
    
    t_ini = tiempo(x_ini);
    t_fin = tiempo(end);
    
    t_interes = tiempo(x_ini:x_fin) - tiempo(x_ini);
    c_interes = salida(x_ini:x_fin);
    
    tiempos_de_interes{end + 1} = t_interes;
    salidas_de_interes{end + 1} = c_interes;

    fprintf("Valor final de respuesta al escalón " + str_amplitud_escalon)
    c_inf = mean(c_interes(end-muestras_para_media, end))
    
    fprintf("Valores de tiempo de establecimiento " + str_amplitud_escalon)
    ks = find(abs(c_interes) >= abs(0.95*c_inf), 1, 'first') - 1
    ts = ks*Ts
    
    figure
    hold on
    stairs(t_interes, abs(c_interes))
    plot([0 t_fin - t_ini], abs(0.95*c_inf)*[1 1], 'r:')
    xlim([0 5])
    title("Respuesta al escalón (" + str_amplitud_escalon + ") de la planta")
    xlabel("Tiempo (s)")
    ylabel("Abs velocidad escalera (m/s) ante escalón de " + str_amplitud_escalon)
    
    fprintf("Modelo de la planta (FdT en m/s/V) (escalon " + str_amplitud_escalon + ")")
    Km = c_inf/amplitud_escalon
    tau = -ts/aproximacion
    
    Kms{end + 1}  = Km;
    taus{end + 1} = tau;

    syms s
    funciones_de_transferencia{end + 1} = Km/(tau*s + 1);

end

fprintf("Funciones de transferencia:")
vpa(funciones_de_transferencia, decimales)

fprintf("Funcion de transferencia media:")
Km  = mean(cell2mat(Kms));
tau = mean(cell2mat(taus));

G = Km/(tau*s + 1);
vpa(G, decimales)

%% 
% 
% 
% *Ejercicio 2*
% 
% Corrobore la validez del modelo mediante su simulación en Matlab/Simulink, 
% como  también se hizo en la Práctica 2
% 
% 

open_system("Fase1/ejercicio2")
set_param("ejercicio2/Escalera mecánica", 'Numerator', string(Km), 'Denominator', '[' + string(tau) + ' 1' + ']')

for i = 1:length(ficheros)
    fichero = ficheros(i).name;
    
    % Saltar '.' y '..'
    if strcmp(fichero, '.') || strcmp(fichero, '..')
        continue;
    end
    
    str_amplitud_escalon = strsplit(fichero, '_'); str_amplitud_escalon = str_amplitud_escalon{end};
    amplitud_escalon = str2double(regexprep(str_amplitud_escalon, '[^-\d]', ''));
    
    set_param("ejercicio2/Consigna (V)", 'After', string(amplitud_escalon))
    sim("ejercicio2")

    t_sim = salidaEscalonMotor(:, 1);
    c_sim = salidaEscalonMotor(:, 2);
    [c_ind, t_ind] = step(amplitud_escalon * tf(Kms{i - 2}, [taus{i - 2} 1]), 5);
    figure
    hold on
    title("Comparación modelo teórico particular, modelo medio del experimento y medidas " + str_amplitud_escalon)
    stairs(t_ind, c_ind)
    stairs(t_sim, c_sim)
    stairs(tiempos_de_interes{i - 2}, salidas_de_interes{i - 2})
    legend('Modelo individual', 'Modelo medio', 'Medidas reales','Location','southeast')
    xlim([0 5])
    xlabel("Tiempo (s)")
    ylabel("Velocidad escalera (m/s)")
end
%% 
% 
% 
% *Ejercicio 3*
% 
% Proponga e implemente los siguientes experimentos para caracterizar las alinealidades 
% de  la planta: 
%% 
% # Para modelar la saturación excite el SM con grandes tensiones 
% # Para modelar  la zona muerta excite el SM con pequeñas tensiones. 
%% 
% Debe de obtener los valores adecuados  con sus unidades de la saturación y 
% zona muerta, explicando el procedimiento de obtención,  y añadir al modelo de 
% Simulink los bloques de las alinealidades caracterizadas  convenientemente colocadas 
% y con coherencia de unidades. 

% Valores obtenidos experimentalmente
fprintf("Saturación de salida (Expresada como módulo máximo de la velocidad de la escalera en m/s). Simétrica.")
c_max_abs = 2.11 % medida lab

fprintf("Zona muerta. (Expresada como la tensión mínima para que se comience a mover la escalera). Asimétrica")
deadzone_in_pos = 0.21 % medida lab
deadzone_out_pos = -0.19 % medida lab

% Compensación de Km
c_inf_deadzone = 0.1517; % Valor final con la zona muerta y la antigua FdT
km_comp = Km * c_inf/c_inf_deadzone;

fprintf("Nueva aproximación de la planta teniendo en cuenta la zona muerta [m/s/V]")
syms s
Gcomp = km_comp/(tau*s + 1);
vpa(Gcomp, decimales)
%% 
% 
% 
% 
% 
% *Ejercicio 4*
% 
% Calcule el equivalente discreto de la planta y caracterice teóricamente su 
% respuesta ante  entrada escalón, tanto en régimen permanente como transitorio. 

BoG = c2d(tf(Km, [tau 1]), Ts) % BoG = (1-z^(-1)) * sum_(Polos G(s)/s) Res[G(s)/(s*(1-exp(s*Ts)*z^(-1))]
[~, zp, ~] = zpkdata(BoG, 'v')

fprintf("Valor final ante escalón unitario")
c_disc_inf = dcgain(BoG) % lim(z->1) BoG = 0.02345 / (1 - 0.8779)

fprintf("Tiempo de establecimiento")
ks_disc = ceil(aproximacion / log(abs(zp))) % Aproximación de log(0.05)/log(abs(zp))
ts_disc = Ts * ks_disc
%% 
% 
% 
% *Ejercicio 5*
% 
% Con el número máximo de personas, el escalón unidad no es suficiente para 
% iniciar el movimiento. La perturbación máxima para el valor positivo es de -0.2m/s, 
% valor con módulo superior a la salida frente al escalón unidad (0.192 m/s)
% 
% 
%% Fase 2
% *Ejercicio 6*
% 
% Calcule la FdeT del regulador con sus correspondientes unidades, mediante 
% la aplicación  detallada del método directo, el resultado debe presentarlo de 
% la forma factorizada,  polinomial y con potencias negativas en z. Presente la 
% ecuación en diferencias del  regulador.

clear;

decimales = 5;

% Planta
Ts=0.035;
Kg = 0.2345;
pg = 0.8779;
aproximacion = -3; %log(0.05)

% Especificaciones deseadas
ts = 0.980;
erpp = 0.03;

fprintf("Parámetros del sistema. Módulo del polo y Km para cada signo")
pm_mod = exp(Ts*aproximacion/ts) 
Km_pos = (1 - pm_mod) * (1 - erpp)
Km_neg = (1 + pm_mod) * (1 - erpp)

fprintf("Funciones de transferencia del controlador")
syms z
BoG = Kg/(z - pg);
F_pos = (Km_pos * (z - pg)) / (Kg * (z - pm_mod - Km_pos));
F1 = vpa(F_pos, decimales)
F_neg = (Km_neg * (z - pg)) / (Kg * (z + pm_mod - Km_neg));
F2 = vpa(F_neg, decimales + 1) 
%% 
% Solución:
% 
% $$F_1(z) = \frac{0.09855z - 0.08652}{0.2345z-0.2338} = \frac{0.4203 - 0.3690\cdot 
% z^{-1}}{1-0.9970\cdot z^{-1}} = 0.4203\cdot\frac{z-0.8779}{z-0.9970}$$
% 
% $$c[k] = 0.9970\cdot r[k-1] + 0.4203\cdot r[k] - 0.3690\cdot r[k-1]$$
% 
% 
% 
% $$F_2(z) = \frac{1.8414z-1.6166}{0.2345z-0.22114}=\frac{7.8525-6.8938\cdot 
% z^{-1}}{1-0.9430\cdot z^{-1}}=7.8525\cdot\frac{z-0.8779}{z-0.9430}$$
% 
% $$c[k] = 0.9430\cdot c[k-1] + 7.8525\cdot r[k] - 6.8938\cdot r[k-1]$$
% 
% *Ejercicio 7*
% 
% Dibuje una topología (diagrama de bloques) del sistema de acuerdo a las unidades 
% elegidas,  que deberán aparecer indicadas en cada señal o secuencia del dibujo. 

fprintf("Funciones de transferencia del sistema con los reguladores")
function res = normalize_to_unit_Z(f, z)
    [num, den] = numden(f);

    coeffs_den = coeffs(den, z);
    z_coeff = coeffs_den(2);
    normalized_num = num / z_coeff;
    normalized_den = den / z_coeff;

    res = vpa(normalized_num / normalized_den, 5);
end

M_pos = normalize_to_unit_Z(simplifyFraction(BoG*F_pos/(1+BoG*F_pos)), z)
M_neg = normalize_to_unit_Z(simplifyFraction(BoG*F_neg/(1+BoG*F_neg)), z)
%% 
% *Ejercicio 8*
% 
% Calcule la FdeT del sistema y caracterice teóricamente su salida ante entradas 
% escalón de  ±1 m/s. 

pm_pos = 0.8984;
km_pos = 0.098555;

pm_neg = -0.8984;
km_neg = 1.8414;

% Queda una muestra de más en ts xq se usa la aproximación y por solo tener
% cuatro decimales en el polo. Si se usan todos los decimales o se utiliza
% log(0.05) se obtiene el resultado exacto
fprintf("Caracterización Mpos")
ts_pos = Ts * ceil(aproximacion/log(abs(pm_pos)))
c_inf_pos = km_pos / (1 - pm_pos)
erpp_pos = 1 - c_inf_pos

fprintf("Caracterización Mneg")
ts_neg = Ts * ceil(aproximacion/log(abs(pm_neg)))
c_inf_neg = km_neg / (1 - pm_neg)
erpp_neg = 1 - c_inf_neg
%% 
% *Ejercicio 9*
% 
% Calcule en régimen permanente tanto la acción de control como la velocidad 
% de la escalera  subiendo y bajando en diferentes condiciones de carga de la 
% escalera, para consignas de  ±1 m/s. Relacione los incrementos o decrementos 
% de acción de control con el peso de las  personas presentes en la escalera. 

pert_5personas = 0.1;
pert_10personas = 0.2;
entrada = 1;
%% 
% Aplicando teorema de superposición

M_r_pos = M_pos;
M_p_pos = normalize_to_unit_Z(-1/(1+BoG*F_pos), z)

M_r_neg = M_neg;
M_p_neg = normalize_to_unit_Z(-1/(1+BoG*F_neg), z)

K0_r_pos = round(subs(M_r_pos, z, 1), decimales)
K0_p_pos = round(subs(M_p_pos, z, 1), decimales)

K0_r_neg = round(subs(M_r_neg, z, 1), decimales)
K0_n_neg = round(subs(M_p_neg, z, 1), decimales)
%% 
% La perturbación no afecta al tiempo de establecimiento, solo al valor final 
% (y por tanto al error)

K0_ref = subs(M_r_pos, z, 1);
K0_pert= subs(M_p_pos, z, 1);

salida_5personas = round(K0_ref * entrada + K0_pert * pert_5personas, decimales)
erpp_5personas = round((entrada - K0_ref * entrada - K0_pert * pert_5personas)/entrada, decimales)

salida_10personas = round(K0_ref * entrada + K0_pert * pert_10personas, decimales)
erpp_10personas = round((entrada - K0_ref * entrada - K0_pert * pert_10personas)/entrada, decimales)

%% Fase 3
% *Ejercicio 10*
% 
%
