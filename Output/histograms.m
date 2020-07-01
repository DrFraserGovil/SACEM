tab = readtable("VeryLong/SuccessGrid.dat");

tab = table2array(tab);




cla;
%colormap(flipud(gray));
image([0,1],[0,20],transpose(tab(2:end,2:end)),'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
c = colorbar;
c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";
%caxis([10^-4 1])
%set(gca,'colorscale','log')