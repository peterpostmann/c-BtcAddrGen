/* Copyright 2014, Kenneth MacKay. Licensed under the BSD 2-clause license. */

#ifndef _MICRO_ECC_H_
#define _MICRO_ECC_H_

#include <stdint.h>

/* uECC_SQUARE_FUNC - If enabled (defined as nonzero), this will cause a specific function to be
used for (scalar) squaring instead of the generic multiplication function. This will make things
faster by about 8% but increases the code size. */
#ifndef uECC_SQUARE_FUNC
    #define uECC_SQUARE_FUNC 1
#endif

#define ECC_KEY_SIZE 32
#define ECC_KEY_PRIVATE_SIZE ECC_KEY_SIZE
#define ECC_KEY_PUBLIC_SIZE (2 * ECC_KEY_SIZE)

typedef struct EccPrivateKey {
    unsigned char key[ECC_KEY_SIZE];
} EccPrivateKey;

typedef struct EccPublicKey {
    unsigned char x[ECC_KEY_SIZE];
    unsigned char y[ECC_KEY_SIZE];
} EccPublicKey;

#ifdef __cplusplus
extern "C"
{
#endif

/* uECC_compute_public_key() function.
Compute the corresponding public key for a private key.

Inputs:
    private_key - The private key to compute the public key for

Outputs:
    public_key - Will be filled in with the corresponding public key

Returns 1 if the key was computed successfully, 0 if an error occurred.
*/
int uECC_compute_public_key(const uint8_t private_key[ECC_KEY_SIZE],
    uint8_t public_key[ECC_KEY_SIZE * 2]);

#ifdef __cplusplus
} /* end of extern "C" */
#endif

#endif /* _MICRO_ECC_H_ */
