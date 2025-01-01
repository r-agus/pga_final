function get_eq_diff(transfer_function, decimales, input, output)
% Funcion para generar la ecuacion en diferencias a partir 
% de la funcion de transferencia. Genera la forma matematica
% y el formato para la GUI.
    [n, d] = tfdata(transfer_function, 'v');
    orden  = length(d);
    str    = "c[k] =";
    str_c  = "y0 = ";
    
    % Normalizar a d(1) unitario
    n = n / d(1);
    d = d / d(1);

    if n(1) > 0
        str = str + sprintf("%.0" + decimales + "f*" + input + "[k]", n(1));
        str_c = str_c + sprintf("%.0" + decimales + "f*x0", n(1));
    else 
        if n(i) < 0
            str = str + sprintf(" - %.0" + decimales + "f*" + input + "[k]", abs(n(1)));
            str = str + sprintf(" - %.0" + decimales + "f*x0", abs(n(1)));
        end
    end
    for i = 2:orden
        if n(i) > 0
            str = str + sprintf(" + %.0" + decimales + "f*" + input + "[k-%d]", n(i), i-1);
            str_c = str_c + sprintf(" + %.0" + decimales + "f*x%d", n(i), i-1);
        else 
            if n(i) < 0
                str = str + sprintf(" - %.0" + decimales + "f*" + input + "[k-%d]", abs(n(i)), i-1);
                str_c = str_c + sprintf(" - %.0" + decimales + "f*x%d", abs(n(i)), i-1);
            end
        end
    end    

    for i = 2:orden
       if d(i) > 0
           str = str + sprintf(" - %.0" + decimales + "f*" + output + "[k-%d]", d(i), i-1);
           str_c = str_c + sprintf(" - %.0" + decimales + "f*y%d", d(i), i-1);
       else 
           if d(i) < 0
                str = str + sprintf(" + %.0" + decimales + "f*" + output + "[k-%d]", abs(d(i)), i-1);
                str_c = str_c + sprintf(" + %.0" + decimales + "f*y%d", abs(d(i)), i-1);
           end
       end
    end

    fprintf("Eq en diferencias: %s\nEq para GUI: %s;", str, str_c)
end
