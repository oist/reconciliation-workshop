#!/bin/bash

echo "# Comparing AleRax and ALE results"

echo "## Comparing logl"
for d in {red,green,blue,rgb}; do
	echo "### Dataset: ${d}"
	#diff ${d}Nov/per_fam_likelihoods.txt ${d}April/per_fam_likelihoods.txt
	#echo "ALE vs AleRax:"
	ALE_LOGL=`grep ">logl" S.tree_${d}*.uml_rec | sed 's/>logl: //'`
	ALERAX_LOGL=`cat ${d}April/per_fam_likelihoods.txt | cut -d " " -f 2`
	printf "| Software | logL |\n|-|-|\n| ALE | ${ALE_LOGL} |\n|AleRax | ${ALERAX_LOGL} |\n"
done


echo "## Comparing DTL rates"
for d in {red,green,blue,rgb}; do
	echo "### Dataset: ${d}"
	#diff ${d}Nov/model_parameters/model_parameters.txt ${d}April/model_parameters/model_parameters.txt
	#echo "ALE vs AleRax:"
	printf "| Software | D | T | L |\n|-|-|-|-|\n|ALE|"
	grep "rate of" -A 1 S.tree_${d}*.uml_rec | tail -n 1 | awk '{print $2"|"$3"|"$4"|"}'
	printf "|AleRax|"
	head -n 2 ${d}April/model_parameters/model_parameters.txt | tail -n 1 | awk '{print $2"|"$4"|"$3"|"}'
done


echo "## Comparing reconciled event numbers"
for d in {red,green,blue,rgb}; do
	echo " "
	echo "### Dataset: ${d}"
	#diff ${d}Nov/reconciliations/perspecies_eventcount.txt ${d}April/reconciliations/totalSpeciesEventCounts.txt 
	#echo "AleRax Nov"
	#cat ${d}Nov/reconciliations/perspecies_eventcount.txt 
	echo "#### AleRax April"
	head -n 1 ${d}April/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/, / /g' | awk '{print "|"$1"|"$2"|"$3"|"$5"|"$4"|"$6"|"$7"|"$8"|"$9"|"}'
	echo "|-|-|-|-|-|-|-|-|-|"
	tail -n+2 ${d}April/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/, / /g' | awk '{print "|"$1"|"$2"|"$3"|"$5"|"$4"|"$6"|"$7"|"$8"|"$9"|"}'
	echo " "
	echo "#### ALE"
## of	 		Duplications	Transfers	Losses	Originations	copies	singletons	extinction_prob	presence	LL
#S_terminal_branch	A(0)	0	0.1915	0.1855	0.0108	1	0.8263	0.110227	1	-13.8975
	printf "|node|duplications|transfers|losses|originations|copies|singletons|presence|\n|-|-|-|-|-|-|-|-|-|\n"
	tail -n11  S.tree_${d}*.uml_rec | sed -E 's/\t/ /' | awk '{print "|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$11"|"}'
	echo " " 
done



echo "## Comparing reconciled event numbers if AleRax is constraied to the ML rates of ALE"
for d in {red,green,blue}; do
	echo "### Dataset: ${d}"

	echo "#### AleRax April"
	head -n 1 ${d}April_fixed/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/, / /g' | awk '{print "|"$1"|"$2"|"$3"|"$5"|"$4"|"$6"|"$7"|"$8"|"$9"|"}'
	echo "|-|-|-|-|-|-|-|-|-|"
	tail -n+2 ${d}April_fixed/reconciliations/totalSpeciesEventCounts.txt | sed -E 's/, / /g' | awk '{print "|"$1"|"$2"|"$3"|"$5"|"$4"|"$6"|"$7"|"$8"|"$9"|"}'
	echo " "
	echo "#### ALE"
	printf "|node|duplications|transfers|losses|originations|copies|singletons|presence|\n|-|-|-|-|-|-|-|-|-|\n"
	tail -n11  S.tree_${d}*.uml_rec | sed -E 's/\t/ /' | awk '{print "|"$2"|"$3"|"$4"|"$5"|"$6"|"$7"|"$8"|"$9"|"$11"|"}'
	echo " "
done

