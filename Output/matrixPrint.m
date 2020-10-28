function matrixPrint(f1,f2,root,divider,mode,sourceRoot)
  
    fracCode = "a";
    if mode(1) == "All"
        fracCode = "b";
    end

    modeName = fracCode + string(mode(2));
   
    fileRoot = sourceRoot + "Model" + root + "_Fraction" + mode(1) + "_Gradient" + mode(2) + "/";
    nModels = readmatrix(fileRoot + "Progress.dat");
    pow = floor(log10(nModels));
    models = "";
    if pow > 3
        base = round(nModels/10^(pow - 1))/10;
        if base > 1
            models = num2str(base) + "$\times";
        else
            models = "$";
        end
        models = models + "10^{" + num2str(pow) + "}$";
    else
        models = num2str(nModels);
    end
    
    %success = readmatrix(root + "SuccessCounts.dat","Delimiter",",");
    tab1 = readmatrix(fileRoot + f1 + ".dat","Delimiter",",");
    tab2 = readmatrix(fileRoot + f2 + ".dat","Delimiter",",");

    if all(size(tab1) == size(tab2)) == 0
        return
    end
    %tab=success;

    if divider == true
        tab = tab1./tab2;
    else
        tab = tab1;
    end
    %         N = readmatrix(root + "Progress");
    %         r = N*length(tab)^2;
    %         q = sum(sum((success)));
    %         successFrac = round(100*q/r,2);
    %         fprintf("%d Galaxies generated for %d models.\n %d Successful models for %f %% success\n", N,r,q,successFrac)
    %cla;

    %colormap(flipud(gray));

    oldMin = min(min(tab));
    mod = 0.95;
    if oldMin < 0
        mod = 1.05;
    end
    newMin = oldMin * mod;
 
    newMax = max(max(tab));
    if newMin == 0
        newMin = -newMax/100;
    end

    I = isnan(tab) | (tab==0);
    tab(I) = newMin - abs(0.05*newMax);
    image([0,1],[0,16],transpose(tab),'CDataMapping','scaled')
    set(gca,'YDir','normal')
    t = "Model \texttt{" + root + modeName+ "}" + " (" + models + " models)";
    title(t,"Interpreter","latex","FontSize",20)
    xlabel("Collapsar Fraction, $f_{coll}$","Interpreter","latex","FontSize",20);
    ylabel("Collapsar Cutoff time, $\tau_{coll}$", "Interpreter","latex","FontSize",20);
    colormap([0 0 0; parula(512)])
    c = colorbar;
    
    
    c.Label.String = f1 ;
    if divider == true
        c.Label.String = f1 + " / " + f2;
    end
    c.Label.FontSize = 17;
    c.Label.Interpreter = "none";

    
    %caxis([newMin, newMax]);

    
    set(gca,'colorscale','linear')

end