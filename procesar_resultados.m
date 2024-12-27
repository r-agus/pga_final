function [ts, error, c_inf] = procesar_resultados(ruta_fichero, Ts, muestras_para_media, decimales, salida_simulacion)
    nombre_fichero = split(ruta_fichero, "/"); nombre_fichero = nombre_fichero(end);

    amplitud_perturbacion = str2num(regexprep(nombre_fichero, "escalon_(-?\d+)m_s_(\d+)_.*", "$1 $2"));
    amplitud = amplitud_perturbacion(1);
    personas = amplitud_perturbacion(2);
    
    contenido = load(ruta_fichero);

    tiempo     = contenido(:, 1);
    referencia = contenido(:, 2);
    salida     = contenido(:, 3);

    x_ini = find(abs(referencia) > 0, 1, 'first');  
    x_fin = length(referencia);

    t_ini = tiempo(x_ini);
    t_fin = tiempo(end);

    t_interes = tiempo(x_ini:x_fin) - tiempo(x_ini);
    c_interes = salida(x_ini:x_fin);
    
    if nargin == 5
        salida_simulacion(end+1:length(t_interes)) = salida_simulacion(end);
    end

    c_inf = mean(c_interes(end-muestras_para_media, end));

    ks = find(abs(c_interes) >= abs(0.95*c_inf), 1, 'first') - 1;

    figure
    hold on
    stairs(t_interes, c_interes)
    if nargin == 5
        stairs(t_interes, salida_simulacion)
    end
    plot([0 t_fin - t_ini], 0.95*c_inf*[1 1], 'r:')
    %xlim([0 5])
    ylim(1.1 * [min(0, c_inf) max(0, c_inf)])
    sgtitle("Respuesta al escalón de amplitud " + amplitud + " m/s con " + personas + " personas de la planta")
    xlabel("Tiempo (s)")
    ylabel("Velocidad escalera (m/s)")

    ts = ks*Ts;
    error = c_inf - amplitud;

    fprintf("Escalón de %d m/s, con %d personas.\n\tMedia de las últimas %d muestras = %.0" + decimales + "f m/s\n\tts  = %.0" ...
        + decimales + "f s\n", amplitud, personas, muestras_para_media, c_inf, ts);
end