#!/bin/bash

echo "Comparing logl"
for d in {red,green,blue,rgb}; do
	echo "dataset: ${d}"
	#diff ${d}Nov/per_fam_likelihoods.txt ${d}April/per_fam_likelihoods.txt
	echo "ALE vs AleRax:"
	grep ">logl" S.tree_${d}*.uml_rec
	cat ${d}April/per_fam_likelihoods.txt
done

echo "======================"

echo "Comparing model parameters"
for d in {red,green,blue,rgb}; do
	echo "dataset: ${d}"
	#diff ${d}Nov/model_parameters/model_parameters.txt ${d}April/model_parameters/model_parameters.txt
	echo "ALE vs AleRax:"
	grep "rate of" -A 1 S.tree_${d}*.uml_rec
	head -n 2 ${d}April/model_parameters/model_parameters.txt
done

echo "======================="

echo "Comparing reconciled event numbers"
for d in {red,green,blue,rgb}; do
	echo " "
	echo "dataset: ${d}"
	#diff ${d}Nov/reconciliations/perspecies_eventcount.txt ${d}April/reconciliations/totalSpeciesEventCounts.txt 
	#echo "AleRax Nov"
	#cat ${d}Nov/reconciliations/perspecies_eventcount.txt 
	echo "AleRax April"
	cat ${d}April/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/ /\t/g'
	echo " "
	echo "ALE"
	tail -n12  S.tree_${d}*.uml_rec | sed -E 's/Duplications/\t\tDuplications/'
	echo "========================="
	echo " " 
done



echo "Comparing reconciled event numbers if AleRax is constraied to the ML rates of ALE"
for d in {red,green,blue}; do
	echo "dataset: ${d}"
	echo "AleRax April"
	cat ${d}April_fixed/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/ /\t/g'
	echo " "
	echo "ALE"
	tail -n12  S.tree_${d}*.uml_rec | sed -E 's/Duplications/\t\tDuplications/'
	echo "========================="
	echo " " 
done

