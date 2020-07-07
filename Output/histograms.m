tab = readtable("Tester/SuccessGrid.dat");

tab = table2array(tab);




cla;
%colormap(flipud(gray));
image([0,1],[0,20],transpose(tab),'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
c = colorbar;
c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";
caxis([0 20])
set(gca,'colorscale','linear')