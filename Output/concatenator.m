root = "StellarParameters/";
subs = "Attempt*";
subdirs = dir(root + subs);
names = {subdirs.name}

A = [];
reshaped = false;
%colormap(flipud(gray));
for name = names
    f = readmatrix(strcat(root + name, "/SuccessGrid.dat"));
    if reshaped == false
        [x,y] = size(f);
        A = zeros(x,y);
        reshaped = true;
    end
    A = A + f;
end
A = transpose(A);

tauMax = 20;

image([0,1],[0,tauMax],A,'CDataMapping','scaled')
set(gca,'YDir','normal')
title("Model Counting","Interpreter","latex","FontSize",20)
xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
ylabel("Collapsar Cutoff time, $\tau_{coll}$ (Gyr)", "Interpreter","latex","FontSize",20);
c = colorbar;
c.Label.String = "Number of Successful Models";
c.Label.FontSize = 24;
c.Label.Interpreter = "Latex";
%caxis([10^-4 1])
set(gca,'colorscale','linear')