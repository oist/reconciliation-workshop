/*
 Copyright (C) 2015-21 Diego Darriba, Alexey Kozlov

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "corax/corax.h"
#include "topological_tests.h"

#define DEBUG_MODE 0

static int kh_core(double **persite_lnl_1, 
                    double **persite_lnl_2,
                    int *partition_lens, 
                    int partition_count, 
                    int nBootstrap,
                    double p_value_threshold)

{   
    int i,j;
    int retval = CORAX_SUCCESS;

    // calculating likelihoods
    double lnl1=0, lnl2=0;

    for (i=0; i<partition_count; i++){
        for (j=0; j<partition_lens[i]; j++){
            lnl1 += persite_lnl_1[i][j];
            lnl2 += persite_lnl_2[i][j];
        }
    }


    // delta Lnl
    double Delta_Lnl = (lnl1 - lnl2);

    if (Delta_Lnl < -1e-3) // to avoid numerical errors
    {
        corax_set_error(CORAX_KH_ORDER_ERROR,
                        "The first array (**persite_lnl_1) should contain the persite LH of the ML tree.\n"
                        "Please change the order of the given arrays, since it turns out that the second tree has higher LH score.\n");
        printf("\nERROR: %s", corax_errmsg);
        return CORAX_FAILURE;
    }

    if (Delta_Lnl < 0) Delta_Lnl = 0;

    // calculating total len
    int nSites = 0;
    for (i = 0; i<partition_count; i++) nSites += partition_lens[i];

    unsigned int _seed = 0;
    corax_random_state* rstate = corax_random_create(_seed);

    // RELL bootstrap process
    double Delta_Ls[nBootstrap];
    double mean_Dl = 0;
    
    for (i = 0; i<nBootstrap; i++){
        
        int pIndex, sIndex;
        lnl1=0;
        lnl2=0;

        for(j = 0; j<nSites; j++){

            if(partition_count > 1){
                pIndex = corax_random_getint(rstate, partition_count);
            } else {
                pIndex = 0;
            }

            sIndex = corax_random_getint(rstate, partition_lens[pIndex]);
            lnl1 += persite_lnl_1[pIndex][sIndex];
            lnl2 += persite_lnl_2[pIndex][sIndex];
        }

        Delta_Ls[i] = lnl1 - lnl2;
        if (DEBUG_MODE) printf("Bootstrap #%d,\tLH1 - LH2 = %f\n",i+1, Delta_Ls[i]);

        mean_Dl += Delta_Ls[i] / nBootstrap;
    }
    
    // substracting the mean and counting the number of samples worse that Delta_Lnl
    int KHsupport=0;
    for(i=0; i<nBootstrap; i++) {
        
        Delta_Ls[i] = Delta_Ls[i] - mean_Dl;
        if ( Delta_Ls[i] >= Delta_Lnl ) KHsupport++;
    }

    double p_value = (double) KHsupport / nBootstrap    ;
    if (DEBUG_MODE) printf("p-value %f\n", p_value);

    if (p_value <= p_value_threshold)
        retval = CORAX_FAILURE;

    return retval;
}


CORAX_EXPORT int corax_KH_RELL_test(double **persite_lnl_1, 
                                    double **persite_lnl_2,
                                    int *partition_lens, 
                                    int partition_count, 
                                    int nBootstrap,
                                    double p_value_threshold)
{
    
    return kh_core(persite_lnl_1,
                    persite_lnl_2,
                    partition_lens,
                    partition_count,
                    nBootstrap,
                    p_value_threshold);
}

