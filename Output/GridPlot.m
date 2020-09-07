



%files = ["b1","b2","BaseCooling","ColdStarRatio","Coll_CoolFrac","Delta_Colls","EuFe_Inf",...,
 %   "f_CCSN","f_Coll","f_NSM","f_SNIa","FeH_Inf","HotStarRatio","M0","M1","M2","MgFe_0","MgFe_Inf","Mu_Stellar","NSM_CoolFrac",...,
 %   "nu_modifier","nu_NSM","nu_SFR","nu_SNIa","OutflowFrac","SFR","s-Fraction","SNIa_CoolFrac","StarMassRatio","tau_NSM","tau_SNIa","X"];
files = ["b1","b2"];
F = figure(1);
sourceRoot = "JointModels/";
 set(F,'units','normalized','outerposition',[0 0 1 1]);
for a = 1 :length(files)-1
    
    f1 = files(a);
    for b = a+1:length(files)
        f2 = files(b);
        f1 + "  " + f2
        tilePrint(f1,f2,true,sourceRoot);
        
        if (a == 1)
           tilePrint(f2,f2,false,sourceRoot);
           if (b==2)
               tilePrint(f1,f1,false,sourceRoot);
               tilePrint("SuccessCounts","SuccessCounts",false,sourceRoot);
           end
        end
        
    end
end



