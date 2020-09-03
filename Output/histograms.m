



files = ["b1","b2","BaseCooling","ColdStarRatio","Coll_CoolFrac","Delta_Colls","EuFe_Inf",...,
    "f_CCSN","f_Coll","f_NSM","f_SNIa","FeH_Inf","HotStarRatio","M0","M1","M2","MgFe_0","MgFe_Inf","Mu_Stellar","NSM_CoolFrac",...,
    "nu_modifier","nu_NSM","nu_SFR","nu_SNIa","OutflowFrac","SFR","s-Fraction","SNIa_CoolFrac","StarMassRatio","tau_NSM","tau_SNIa","X"];
 F = figure('units','normalized','outerposition',[0 0 1 1]);
for a = 1 :length(files)-1
    
    f1 = files(a);
    for b = a+1:length(files)
        f2 = files(b);
        f1 + "  " + f2
        tilePrint(f1,f2,true);
        
        if (a == 1)
           tilePrint(f2,f2,false);
           if (b==2)
               tilePrint(f1,f1,false);
               tilePrint("SuccessCounts","SuccessCounts",false);
           end
        end
        
    end
end


function tilePrint(f1,f2,duoMode)
    cla;
    roots = ["L", "M", "T", "S"];
    modes = ["All","Cold"];
    for j = 1:length(modes)
        mode = modes(j);

       

        for i = 1:length(roots)
            subplot(2,2,i);
            matrixPrint(f1,f2,roots(i),duoMode,mode);

            

        end
          svT = "Images/"+mode+"/";
         
        if duoMode == true
            svT = svT + f1 + "_AND_" + f2;
        else
            svT = svT + "just_" + f1;
        end
        
        svT = svT + ".png";
        
        saveas(gca,svT);
        

    end
end
function matrixPrint(f1,f2,root,divider,mode)

    modeName = "a";
    if mode == "Cold"
        modeName = "b";
    end

    fileRoot = "Model" + root + "_" + mode + "/";
    %success = readmatrix(root + "SuccessCounts.dat","Delimiter",",");
    tab1 = readmatrix(fileRoot + f1 + ".dat","Delimiter",",");
    tab2 = readmatrix(fileRoot + f2 + ".dat","Delimiter",",");

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
    t = "Model \texttt{" + root + modeName+ "}";
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
