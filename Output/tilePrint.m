function tilePrint(f1,f2,duoMode,sourceRoot)
    cla;
    roots = ["L", "M", "T", "S"];
    fracmodes = ["All","Cold"];
    gradmodes = ["0","1","2","3"];
    for j = 1:length(fracmodes)
        for k = 1:length(gradmodes)
            mode = [fracmodes(j), gradmodes(k)];



            for i = 1:length(roots)
                subplot(2,2,i);
                matrixPrint(f1,f2,roots(i),duoMode,mode,sourceRoot);



            end
            
            fracCodes = ["b","a"];
            rooter = fracCodes(j) + gradmodes(k);
            svT = "Images/Mode_"+rooter+"/";
            if ~exist(svT,'dir')
                mkdir(svT);
            end

            
            if duoMode == true
                reroot1 = svT + f1 + "Grids/";
                reroot2 = svT + f2 + "Grids/";
                fName1 = f1 + "_AND_" + f2;
                fName2 = f2 + "_AND_" + f1;
            else
                reroot1 = svT + f1 + "Grids/";
                reroot2 = svT + "SingleVariableGrids/";
                fName1 = "just_" + f1;
                fName2 = f1;
                
                if f1 == "SuccessCounts"
                    reroot2 = "Images/ExclusionGrids/";
                    fName2 = "model_" + rooter;
                end
            end
            
            rerooters = [reroot1, reroot2];
            fNames = [fName1, fName2];
            
            for i = 1:length(rerooters)
            
                if ~exist(rerooters(i),'dir')
                    mkdir(rerooters(i));
                end
                svT = rerooters(i) + fNames(i) + ".png";
                saveas(gca,svT);
            end
        end
    end
end
