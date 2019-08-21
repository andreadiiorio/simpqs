/*
 * gmp_patch.h
 *
 *  Created on: Dec 26, 2011
 *      Author: martani
 */

#include <gmp.h>

#ifndef GMP_PATCH_H_
#define GMP_PATCH_H_


#define LEGANDRE_EXTRA_CHECK_SQRTM
int mpz_sqrtm(mpz_t q, const mpz_t n, const mpz_t p);   //TODO SHANKS_TONNARELLI

#endif /* GMP_PATCH_H_ */


