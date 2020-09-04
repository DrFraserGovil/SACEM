#!/bin/bash
for model in L M T S
do
	i=0
	for fractionMode in Cold All
	do
		for gradientMode in 0 1 2 3
		do 
			base="./submitShort -threads 16"
			if [ $gradientMode -eq 3 ]
			then
				base="./submitLong -threads 28"
			fi
			command="${base} -random 10000000 -grid 101 -dir Output/Model${model}_Fraction${fractionMode}_Gradient${gradientMode} -mode 1 -dt 0.2 -constrain ${model} -allfrac ${i} -gradient ${gradientMode}"
			echo $command

		done
		i=$((i + 1))
		echo $i
	done
done
 #./submitLong -threads 28 -random 10000000 -grid 101 -dir Output/ModelS_Cold_Grad -mode 1 -dt 0.2 -constrain s -allfrac 0 -gradient 1
