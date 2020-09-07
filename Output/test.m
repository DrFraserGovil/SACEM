M0 = 8.5;
M1 = 4;
M2 = 25;
b1 = 0.3;
b2 = 14;
nuSFR = 0.5;
t = [0:0.01:40];

base = @(tp,M,b)  M/(nuSFR*b - 1) * (exp(-tp/b) - exp(-nuSFR * tp));



Mc = M0 * exp(-nuSFR*t) + base(t,M1,b1) + base(t,M2,b2);

plot(t,Mc)