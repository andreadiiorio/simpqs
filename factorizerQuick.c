//
// Created by andysnake on 22/08/19.
//

#include <stdlib.h>
#include "factorizerQuick.h"

int factorizeTrialDivide(struct ArrayEntry *arrayEntry) {
    /*
     * trial divide array entry concurrently
     * job enqueue --->
     *  TODO BLOCK ENQUEUE TO AMMORTIZE LOCK/UNLOCK OVERHEAD
     *              <-- job dequeue (factorizerGroupManager)
     *                      concurrent factorize trigger for job
     *                          -->factorizer thread worker new N notify ---> thread_signal vs polling vs condvar
                                  result write inplace arrayEntry

     */
    return EXIT_SUCCESS;
}
