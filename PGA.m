%% Proyecto Global de Asignatura
%% Rubén A. y David A. - Otoño 2024

clear;
decimales = 5;
muestras_para_media = 50;
aproximacion = -3;
%% Fase 1 - Caracterización de la planta
% *Ejercicio 1*
% 
% Capture y caracterice la salida de la planta en función de la tensión de entrada 
% con un procedimiento análogo al empleado en la Práctica 2 y obtenga su función 
% de transferencia (FdeT) en las unidades que considere más adecuadas para trabajar 
% en el desarrollo de todo el diseño. 

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
% *Ejercicio 2*
% 
% Corrobore la validez del modelo mediante su simulación en Matlab/Simulink, 
% como también se hizo en la Práctica 2

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
% *Ejercicio 3*
% 
% Proponga e implemente los siguientes experimentos para caracterizar las alinealidades 
% de la planta: 
%% 
% # Para modelar la saturación excite el SM con grandes tensiones 
% # Para modelar la zona muerta excite el SM con pequeñas tensiones. 
%% 
% Debe de obtener los valores adecuados con sus unidades de la saturación y 
% zona muerta, explicando el procedimiento de obtención, y añadir al modelo de 
% Simulink los bloques de las alinealidades caracterizadas convenientemente colocadas 
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
% *Ejercicio 4*
% 
% Calcule el equivalente discreto de la planta y caracterice teóricamente su 
% respuesta ante entrada escalón, tanto en régimen permanente como transitorio. 

% Se añade un retraso para tener en cuenta el que se quitó en el ejercicio 1
BoG = series(c2d(tf(Km, [tau 1]), Ts), tf(1, [1 0], Ts)) % BoG = ((1-z^(-1)) * sum_(Polos G(s)/s) Res[G(s)/(s*(1-exp(s*Ts)*z^(-1))]) * (1/z)
[~, zp, ~] = zpkdata(BoG, 'v');

fprintf("Valor final ante escalón unitario")
c_disc_inf = dcgain(BoG) % lim(z->1) BoG = 0.02345 / (1 - 0.8779)

fprintf("Tiempo de establecimiento")
ks_disc = ceil(aproximacion / log(abs(max(zp)))) + (length(zp) - 1) % Teniendo en cuenta el retardo
ts_disc = Ts * ks_disc
%% 
% *Ejercicio 5*

pert10personas = 0.2;
%% 
% Con el número máximo de personas, el escalón unidad no es suficiente para 
% iniciar el movimiento positivo. La perturbación máxima para el valor positivo 
% es de -0.2m/s, valor con módulo superior a la salida frente al escalón unidad 
% (0.1993 m/s).

c_inf_pert_pos = max([0, Km*1 - pert10personas])
%% 
% Sin embargo, en el movimiento negativo, la perturbación ""ayuda"" a bajar, 
% por lo que el valor final será de:

c_inf_pert_neg = min([0, Km*(-1) - pert10personas])
%% Fase 2
% *Ejercicio 6*
% 
% Calcule la FdeT del regulador con sus correspondientes unidades, mediante 
% la aplicación detallada del método directo, el resultado debe presentarlo de 
% la forma factorizada, polinomial y con potencias negativas en z. Presente la 
% ecuación en diferencias del regulador.

% Normalizar a planta de la forma Kg/(s-pg) = (Km/tau)/(s - (-1/tau)) = Km/(tau*s + 1)
Kg = Km/tau;
pg = -1/tau;

clearvars -except Kg pg Ts aproximacion decimales;

% Aproximación del BoG, teniendo en cuenta el retardo adicional
%% 
fprintf("Equivalente discreto de la planta a controlar")
G = tf(Kg, [1 -pg]);
retardo = tf(1, [1 0], Ts);
BoG = series(c2d(G, Ts), retardo)
[NumBoG, DenBoG] = tfdata(BoG, 'v')

% Especificaciones deseadas
ts = 0.980;
erpp = 0.03;

fprintf("===== MÉTODO DE TRUXAL =====\n1. Aplicabilidad del modelo. Aplicable por entrada al escalón y polo menor que la unidad\n2. Causalidad")
grado_Den_BoG = length(DenBoG) - find(DenBoG ~= 0, 1, 'first');
grado_Num_BoG = length(NumBoG) - find(NumBoG ~= 0, 1, 'first');
dif_grados_G = grado_Den_BoG - grado_Num_BoG
retardos = dif_grados_G - 1;

fprintf("3. Modelo genérico. Se genera un modelo sobreamortiguado teniendo en cuenta la diferencia de grados del paso anterior\n" ...
    + "\tKm/(z^(-" + retardos + ")*(z-pm))\nCálculo de parámetros:");

ts_prima = ts - retardos*Ts;
pm = exp(Ts * aproximacion / ts_prima) % Ignorando solución negativa por "no ser práctica"
km = (1-pm)*(1-erpp)
M_obj = tf(km, [1, -pm], Ts)*(tf(1, [1 0], Ts)^retardos)

fprintf("4. Cálculo del controlador")
F = zpk(minreal(M_obj/(BoG * (1 - M_obj))))

%% 
% Solución:
% 
% $$F_1(z) = \frac{0.09855z - 0.08652}{0.2345z-0.2338} = \frac{0.4203 - 0.3690\cdot 
% z^{-1}}{1-0.9970\cdot z^{-1}} = 0.4203\cdot\frac{z-0.8779}{z-0.9970}$$
% 
% $$c[k] = 0.9970\cdot r[k-1] + 0.4203\cdot r[k] - 0.3690\cdot r[k-1]$$
% 
% $$F_2(z) = \frac{1.8414z-1.6166}{0.2345z-0.22114}=\frac{7.8525-6.8938\cdot 
% z^{-1}}{1-0.9430\cdot z^{-1}}=7.8525\cdot\frac{z-0.8779}{z-0.9430}$$
% 
% $$c[k] = 0.9430\cdot c[k-1] + 7.8525\cdot r[k] - 6.8938\cdot r[k-1]$$
% 
% *Ejercicio 7*
% 
% Dibuje una topología (diagrama de bloques) del sistema de acuerdo a las unidades 
% elegidas, que deberán aparecer indicadas en cada señal o secuencia del dibujo. 


%% 
% *Ejercicio 8*
% 
% Calcule la FdeT del sistema y caracterice teóricamente su salida ante entradas 
% escalón de ±1 m/s. 

fprintf("Valor de FdeT del sistema controlado")
M_final = zpk(minreal(feedback(series(F, BoG), 1)))
[~, p_m_final, ~] = zpkdata(M_final, 'v');

fprintf("Caracterización de salida frente a escalón unidad positivo")
ts_M = Ts * (ceil(round(aproximacion / log(max(p_m_final)), 10)) + (length(p_m_final) - 1)) % ceil(round(...)) para evitar que x + 1e-12 suba a x + 1
Km_final = dcgain(M_final) % lim_(z->1) M_final
c_inf_mas_1V   = Km_final * 1
c_inf_menos_1V = Km_final * (-1)

%% 
% *Ejercicio 9*
% 
% Calcule en régimen permanente tanto la acción de control como la velocidad 
% de la escalera subiendo y bajando en diferentes condiciones de carga de la escalera, 
% para consignas de ±1 m/s. Relacione los incrementos o decrementos de acción 
% de control con el peso de las personas presentes en la escalera. 

pert_5personas = -0.1;
pert_10personas = -0.2;
%% 
% Aplicando teorema de superposición

M_r = M_final;
M_p = zpk(minreal(feedback(1, series(F, BoG))));

%% 
% La perturbación no afecta al tiempo de establecimiento, solo al valor final 
% (y por tanto al error)

K0_ref = dcgain(M_r);
K0_pert= dcgain(M_p);
K0_F   = dcgain(F);

entradas       = [1 -1];
perturbaciones = [0 -0.1 -0.2];
perturbaciones_names = [" 0 personas", " 5 personas", "10 personas"];

fprintf("Cálculo de error y acción en régimen permanente para distintas perturbaciones\n")

for i = 1:length(entradas)
    fprintf("\nPara escalón de %0.2fV:\n", entradas(i))
    for j = 1:length(perturbaciones)
        r = entradas(i);
        p = perturbaciones(j);

        c = K0_ref * r + K0_pert * p;
        e = c - r;
        a = K0_F * e;
        
        fprintf("\tCon %s:\tc_inf=%.03f\terpp=%.03f\ta_inf=%.03f\n", perturbaciones_names(j), c, e, a);
    end
end

%% Fase 3
% *Ejercicio 10*