function overPlotter(files)
f = readtable("FilteredSAGA.dat");


fe = f{:,1};
mg = f{:,2};
eu = f{:,3};

[fe,I] = sort(fe);
eu = eu(I);
mg = mg(I);
N = ceil(length(fe)/10);
femin = -2.2;


insertCoord = 1;
yEsum = 0;
yMsum = 0;
xsum = 0;
n = 0;
NinBin = N;
for i = 1:length(fe)
	x = fe(i);
	
	yEsum(insertCoord) = yEsum(insertCoord) + eu(i);
	yMsum(insertCoord) = yMsum(insertCoord) + mg(i) - x;
	xsum(insertCoord) = xsum(insertCoord) + x;
	n(insertCoord) = n(insertCoord) + 1;
	if n(insertCoord) >= NinBin
		insertCoord = insertCoord + 1;
		n(end+1) = 0;
		yEsum(end+1) = 0;
		yMsum(end+1) = 0;
		xsum(end+1) = 0;
	end
end
yEBins = yEsum./n;
yMBins = yMsum./n;
xBins = xsum./n;

subplot(2,1,1)

cla;

mgfe = mg - fe;
scatterCut =mgfe > 0.7 | mgfe < -0.2;

scatter(fe(~scatterCut),mgfe(~scatterCut),60,'.')
hold on;
plot(xBins,yMBins,'r-s','MarkerSize',6,'MarkerFaceColor','r')

for i= files
	fileName = i + "/SingleEvaluation.dat";
	f2 = readtable(fileName);
	feFake = f2.Fe_H;
	mgFake = f2.Mg_H;
	euFake = f2.Eu_H;
    plot(feFake,mgFake - feFake);
end
xlabel('[Fe/H]','Interpreter','latex','FontSize',15);
ylabel('[Mg/Fe]','Interpreter','latex','FontSize',15);
%legend({'SAGA database','Bin-Average','Model, $Z_{cut} = 0.1 Z_\odot$','Model, $Z_{cut} = 0.3 Z_\odot$','Model, $Z_{cut} = 0.7 Z_\odot$','Model, $Z_{cut} = Z_\odot$','Model, $Z_{cut} \to \infty$'},'Interpreter','latex','FontSize',15)
%legend('SAGA database','Bin-Average','FontSize',15)
axis([femin,0.5,-0.2,0.7])

subplot(2,1,2)
cla;

scatterCut =eu > 1.3 | eu < -0.2;

scatter(fe(~scatterCut),eu(~scatterCut),60,'.')
hold on;
plot(xBins,yEBins,'r-s','MarkerSize',6,'MarkerFaceColor','r')

for i= files
	fileName = i + "/SingleEvaluation.dat";
	f2 = readtable(fileName);
	feFake = f2.Fe_H;
	mgFake = f2.Mg_H;
	euFake = f2.Eu_H;
    plot(feFake,euFake - feFake);
end
xlabel('[Fe/H]','Interpreter','latex','FontSize',15);
ylabel('[Eu/Fe]','Interpreter','latex','FontSize',15);
legend({'SAGA database','Bin-Average','$Z_{cut} = 0.1 Z_\odot$','$Z_{cut} = 0.3 Z_\odot$', '$Z_{cut} = 0.5 Z_\odot$','$Z_{cut} = 0.7 Z_\odot$','$Z_{cut} = Z_\odot$','$Z_{cut} \to \infty$'},'Interpreter','latex','FontSize',15)
%legend('SAGA database','Bin-Average','SACEM Prediction','FontSize',15)
xlim([femin,0.5])
end
