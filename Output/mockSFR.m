M1 = 8.5;
M2 = 25;
b1 = 0.3;
b2 = 14;
M0 = 10;

nuSFR = 0.005;
nuCool = 0.4;
fc = 0	;


A = @(t) M1/b1*exp(-t/b1) + M2/b2*exp(-t/b2);
SigmaS = 0;
SigmaH = 0;
SigmaC = M0;
SFR = nuSFR*SigmaC;
ts = [0:0.01:10];
i = 1;
for t = ts(1:end-1)
	SigmaS(i+1) = SigmaS(i) +  nuSFR*SigmaC(i);
	SigmaC(i+1) = SigmaC(i) - nuSFR*SigmaC(i) + A(t) + nuCool*SigmaH(i);
	SigmaH(i+1) = SigmaH(i)-nuCool*SigmaH(i);
	SFR(i+1) = SigmaC(i+1);
	i= i+1;
end

subplot(2,1,1)
plot(ts,SFR)

subplot(2,1,2)
cla;
hold on;
plot(ts,SigmaS)
plot(ts,SigmaH)
plot(ts,SigmaC)


Mass = @(t) (10./t).^(2/5);
ZetaT = @(t) IMF(Mass(t));
Integrand = @(t,Q,cg,dt,s) IMF( (10./(t-Q)).^(2/5) ).*Sigma(Q,cg,dt,s)*25.*( (t - Q)/10).^(7/2);


function imf = IMF(M)
	alpha = 2.3;
	
	f = (alpha - 1)/alpha;
	less = (M <= 1)*f;
	more = (M > 1)*f.*(M).^(-alpha);
	
	imf = less + more;
	imf(1) = less(1);
end

function s = Sigma(ts,coldgas,deltaT,start)

	
	s = [];
	for i = 1:length(ts) 
		t = ts(i);
		I = floor((t-start+10*eps)/deltaT);
		root = deltaT*I+start;
		frac = (t - root)/deltaT;
		if abs(frac)<10e-5
			s(i) = coldgas(I+1);
		else

			sD = coldgas(I+1);
			sU = coldgas(I+2);

			s(i) = (sU-sD)*frac + sD;
		end
	end
end

function val = D(t,coldGas)

	ts =   
	p = Integrand(t,ts,cg,cg(1),cg(2)-cg(1));

	trapz(ts,p);
end