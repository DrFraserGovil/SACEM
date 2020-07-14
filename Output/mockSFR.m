M1 = 8;
M2 = 46;
b1 = 0.8;
b2 = 14;
M0 = 0.8;

Masses = [M1,M2];
betas = [1.0/b1, 1.0/b2];

nu = 1.8;
lambda = 1.4;
fh = 0.75	;
mu = 0.03;
dOut = 0.4;

dt = 0.01;
ts = 0:dt:14;

Mc = M0*(1-fh);
Mh = 0;
Ms = M0*fh;

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
plot(ts,Mc*nu./(Mc+Ms));
% plot(ts, Mh);
% plot(ts,Ms);
%plot(ts,Mc+Ms)

f = readtable("massotime.dat");
g = readtable("SFmstat.dat");
nonIn = (f.Var3~=0);
f(nonIn,:) = [];

rT = [];
rC = [];
rH = [];
rS = [];
sfr = [];
for i = 1:height(f)
   t = f.Var1(i);
   if (isempty(rT) || t ~= rT(end) )
      rT(end+1) = f.Var1(i);
      rS(end+1) = f.Var4(i);
      rC(end+1) = f.Var5(i);
      rH(end+1) = f.Var6(i);
      sfr(end+1) = g.Var4(i);
   else
      rS(end) = rS(end) + f.Var4(i);
      rC(end) = rC(end) + f.Var5(i);
      rH(end) = rH(end) + f.Var6(i);
      sfr(end) = sfr(end) + g.Var4(i)*g.Var6(i)/0.015;
   end
   
end
rT = rT*15e-3;
rC = rC/10^10;
rH = rH/10^10;
rS = rS/10^10;
sfr = sfr/10^10;
clear f
clear g
plot(rT,rS)
plot(rT,rC)
%plot(rT,rH)
%plot(rT,rS + rC + rH)
plot(rT,sfr./(rC + rH + rS))
legend({"C","H","S","Total","Total 2"})
ylabel("SFR (10^{10} / Gyr")
xlabel("Time (GYr)")
set(gca,"yscale","linear")