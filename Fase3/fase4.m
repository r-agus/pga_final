clear
aproximacion = -3;
Ts = 0.035;

BoG = zpk(tf(0.02596, [1 -0.8698 0], 0.035))

ts = 1.12;

mod_p = exp(Ts*aproximacion /ts);

figure
rlocus(BoG)

X = minreal(BoG * tf([1 -0.8698], [1 -1], 0.035));
[z,p,k] = zpkdata(X, 'v')

kp = prod(abs(mod_p - p))/(k*prod(abs(mod_p - p)))
%%
figure
hold on
rlocus()
UnitCircle(mod_p)
xlim(1.2*[-1 1])
ylim(1.2*[-1 1])
%%

k = pi; % lo siento mucho (esto es 3.14)
q0 = 1 * k;       % q0 =  kp + ki*Ts/2
q1 = -0.8698 * k; % q1 = -kp + ki*Ts/2

ki = (q0 + q1)/Ts
kp = (q0 - q1)/2

F = tf(kp, 1, Ts) + ki*Ts/2*tf([1 1], [1 -1], Ts)

M = zpk(minreal(feedback(F*BoG, 1)))

[~, p, ~] = zpkdata(M, 'v');
ts = Ts * ceil(aproximacion/log(max(abs(p))))