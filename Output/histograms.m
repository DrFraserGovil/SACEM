tab = readtable("LargeRun.dat");

tab = table2array(tab);
taus = tab(:,1);
matrix = tab(:,2:end);


cla;
%colormap(flipud(gray));
image([0,1],[0,taus(end)],matrix,'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
c = colorbar;
c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";
%caxis([10^-4 1])
set(gca,'colorscale','log')