#!/bin/bash
for model in U W V M
do
	i=0
	for fractionMode in Cold
	do
		for gradientMode in 0 1
		do 
			base="./submitShort -threads 16"
			
			command="${base} -random 10000000 -grid 101 -dir Output_NEW/Model${model}_Fraction${fractionMode}_Gradient${gradientMode} -mode 1 -dt 0.2 -constrain ${model} -allfrac ${i} -gradient ${gradientMode}"
			echo $command

		done
		i=$((i + 1))
		
	done
done
 #./submitLong -threads 28 -random 10000000 -grid 101 -dir Output/ModelS_Cold_Grad -mode 1 -dt 0.2 -constrain s -allfrac 0 -gradient 1
