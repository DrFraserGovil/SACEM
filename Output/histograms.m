
root = "ModelT_All/";
success = readmatrix(root + "SuccessCounts.dat","Delimiter",",");
tab1 = readmatrix(root + "tau_SNIa.dat","Delimiter",",");
tab2 = readmatrix(root + "tau_NSM","Delimiter",",");

tab=success;
tab = tab1;
%tab = tab1-tab2;

N = readmatrix(root + "Progress");
r = N*length(tab)^2;
q = sum(sum(~isnan(success)));
successFrac = round(100*q/r,2);
fprintf("%d Galaxies generated for %d models.\n %d Successful models for %f %% success\n", N,r,q,successFrac)
cla;

%colormap(flipud(gray));
image([0,1],[0,16],transpose(tab),'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
colormap([1 1 1; parula(512)])
c = colorbar;
%c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";

oldMin = min(min(tab));
mod = 0.8;
if oldMin < 0
    mod = 1.05;
end
newMin = oldMin * mod;

newMax = max(max(tab));
if newMin == 0
    newMin = -newMax/100;
end
caxis([newMin, newMax]);
% % 
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