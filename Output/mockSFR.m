M1 = 5;
M2 = 46;
b1 = 0.3;
b2 = 14;
M0 = 0.8;

Masses = [M1,M2];
betas = [1.0/b1, 1.0/b2];

nu = 1;
lambda = 1;
fh = 0.75	;
mu = 0.1;
dOut = 2.5;

dt = 0.0001;
ts = 0:dt:14;

Mc = M0;
Mh = 0;
Ms = 0;

for i = 2:length(ts)
   t = ts(i-1);
   
   
   dAcc = sum(Masses.*betas.*exp(-betas*t));
   dMc = dAcc- (1+dOut)*nu*Mc(i-1) + lambda*Mh(i-1) + (1-fh)*mu*Ms(i-1);
   dMh = fh*mu*Ms(i-1) - lambda*Mh(i-1) + Mc(i-1)*nu*dOut;
   dMs = nu*Mc(i-1) - mu*Ms(i-1);
   
   Mc(i) = Mc(i-1) + dMc * dt;
   Mh(i) = Mh(i-1) + dMh * dt;
   Ms(i) = Ms(i-1) + dMs * dt;
end

cla;
hold on;
plot(ts,Mc);
plot(ts, Mh);
plot(ts,Ms);


r = detectImportOptions("LunchTest/SingleEvaluation.dat");
f = readtable("LunchTest/SingleEvaluation.dat",r);
plot(f.Time,f.Mcg);
plot(f.Time,f.Mhg);
plot(f.Time,f.Ms);


p = mu + nu + lambda*(1+dOut);
q = mu*nu*(fh+dOut)+ lambda*(mu + nu);
zi = -Masses.*(betas.^2 - (mu + lambda)*betas + mu*lambda);
D= mu*lambda*(M0 + sum(Masses))/q;
C = zi./(betas.^2 - p*betas + q);
omega = sqrt(p^2/4 - q);
J = M0 - D - sum(C);
K = sum(betas.*(Masses + C)) - nu*(1+dOut)*M0;

A = J/2 - (p*J + 2*K)/(4*omega);
B = J/2 + (p*J + 2*K)/(4*omega);

mcg = D + exp(-p*ts/2).*(A*exp(-omega*ts) + B*exp(omega*ts) ) + C(1)*exp(-ts/b1) + C(2)*exp(-ts/b2);

plot(ts,mcg+0.01)
legend(["Cold Gas", "Hot Gas", "Stars", "Cold Gas Sim", "Hot Gas Sim", "Stars Sm","Reanalysis"]);