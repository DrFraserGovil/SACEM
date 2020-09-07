function tilePrint(f1,f2,duoMode,sourceRoot)
    cla;
    roots = ["L", "M", "T", "S"];
    fracmodes = ["All","Cold"];
    gradmodes = ["0","1","2"];
    for j = 1:length(fracmodes)
        for k = 1:length(gradmodes)
            mode = [fracmodes(j), gradmodes(k)];



            for i = 1:length(roots)
                subplot(2,2,i);
                matrixPrint(f1,f2,roots(i),duoMode,mode,sourceRoot);



            end
            
            fracCodes = ["a","b"];
            rooter = fracCodes(j) + gradmodes(k);
            svT = "Images/Mode_"+rooter+"/";
            if ~exist(svT,'dir')
                mkdir(svT);
            end
            if duoMode == true
                svT = svT + f1 + "_AND_" + f2;
            else
                svT = svT + "just_" + f1;
            end

            svT = svT + ".png";

            saveas(gca,svT);
        
        end
    end
end
