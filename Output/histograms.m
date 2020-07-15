

tab1 = readmatrix("NewTest/X.dat","Delimiter",",");
tab2 =  readmatrix("NewTest/M2.dat","Delimiter",",");

tab=tab1;
%tab = tab1./tab2;

N = 40000;
r = N*length(tab)^2;
q = sum(sum(~isnan(tab)));
successFrac = round(100*q/r,2);
fprintf("%f %% success\n", successFrac)
cla;
%colormap(flipud(gray));
image([0,1],[0,15],transpose(tab),'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
c = colorbar;
%c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";
%caxis([0 1.5])
% 
% x = 0:0.01:1;
% m = 12;
% inter = 0.45;
% y = m*(x - inter);
% hold on;
% plot(x,y,'r')
% 
% m2 = 11.5;
% inter2 = 0.3;
% y2 = m2*(x -inter2);
% plot(x,y2,'k')
set(gca,'colorscale','linear')