
for model in L M T S
do
	for fractionMode in 0 1
	do
		for gradientMode in 0 1
		do 
			base="./submitLong -threads 28"
			if [ $gradientMode -eq 0 ]
			then
				base="./submitShort -threads 16"
			fi
			command="${base} -random 10000000 -grid 101 -dir Output/Model${model}_${fractionMode}_${gradientMode} -mode 1 -dt 0.2 -constrain ${model} -allfrac ${fractionMode} -gradient ${gradientMode}"
			echo $command
		done
		
	done
done
 #./submitLong -threads 28 -random 10000000 -grid 101 -dir Output/ModelS_Cold_Grad -mode 1 -dt 0.2 -constrain s -allfrac 0 -gradient 1
