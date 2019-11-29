tab = readtable("Data2.dat");
% taus = [];
% fs = [];
% matrix = [];
% N = height(tab);
% wait = 5;
% current = -5;
% for i = 1:N
% 	p = i/N*100;
% 	if p > current+wait
% 		disp(strcat(num2str(round(p,1)),"%"));
% 		current = p;
% 	end
% 	t = tab{i,1};
% 	f = tab{i,2};
% 	
% 	tIndex = find(taus==t);
% 	fIndex = find(fs == f);
% 	newPoint = false;
% 	
% 	if isempty(tIndex)
% 		taus(end+1) = t;
% 		tIndex = length(taus);
% 		newPoint = true;
% 	end
% 	
% 	if isempty(fIndex)
% 		fs(end+1) = f;
% 		fIndex = length(fs);
% 		newPoint = true;
% 	end
% 	
% 	if ~newPoint
% 		matrix(tIndex,fIndex) = matrix(tIndex,fIndex)+1;
% 	else
% 		matrix(tIndex,fIndex) = 1;
% 	end
% 	
% end
tab = table2array(tab);
taus = tab(:,1);
matrix = tab(:,2:end);


cla;
colormap(flipud(gray));
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