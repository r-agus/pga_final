function [str] = get_eq_diff(transfer_function, decimales, input, output)
    [n, d] = tfdata(transfer_function, 'v');
    orden  = length(d);
    str    = "c[k] =";
    
    % Normalizar a d(1) unitario
    n = n / d(1);
    d = d / d(1);
    
    for i = 2:orden
       if d(i) > 0
           str = str + sprintf(" - %.0" + decimales + "f*" + output + "[k-%d]", d(i), i-1);
       else 
           if d(i) < 0
                str = str + sprintf(" + %.0" + decimales + "f*" + output + "[k-%d]", abs(d(i)), i-1);
           end
       end
    end

    if n(1) > 0
        str = str + sprintf(" + %.0" + decimales + "f*" + input + "[k]", n(1));
    else 
        if n(i) < 0
            str = str + sprintf(" - %.0" + decimales + "f*" + input + "[k]", abs(n(1)));
        end
    end
    for i = 2:orden
        if n(i) > 0
            str = str + sprintf(" + %.0" + decimales + "f*" + input + "[k-%d]", n(i), i-1);
        else 
            if n(i) < 0
                str = str + sprintf(" - %.0" + decimales + "f*" + input + "[k-%d]", abs(n(i)), i-1);
            end
        end
    end
end