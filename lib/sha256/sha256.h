/* sha256.h - TinyCrypt interface to a SHA-256 implementation */

/*
 *  Copyright (C) 2017 by Intel Corporation, All Rights Reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *    - Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *
 *    - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 *    - Neither the name of Intel Corporation nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @file
 * @brief Interface to a SHA-256 implementation.
 *
 *  Overview:   SHA-256 is a NIST approved cryptographic hashing algorithm
 *              specified in FIPS 180. A hash algorithm maps data of arbitrary
 *              size to data of fixed length.
 *
 *  Security:   SHA-256 provides 128 bits of security against collision attacks
 *              and 256 bits of security against pre-image attacks. SHA-256 does
 *              NOT behave like a random oracle, but it can be used as one if
 *              the string being hashed is prefix-free encoded before hashing.
 *
 *  Usage:      1) call tc_sha256_init to initialize a struct
 *              tc_sha256_state_struct before hashing a new string.
 *
 *              2) call tc_sha256_update to hash the next string segment;
 *              tc_sha256_update can be called as many times as needed to hash
 *              all of the segments of a string; the order is important.
 *
 *              3) call tc_sha256_final to out put the digest from a hashing
 *              operation.
 */

#ifndef __TC_SHA256_H__
#define __TC_SHA256_H__

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

#ifndef NULL
#define NULL ((void *)0)
#endif

#define TC_CRYPTO_SUCCESS 1
#define TC_CRYPTO_FAIL 0

#define TC_ZERO_BYTE 0x00

/**
* @brief Copy the the buffer 'from' to the buffer 'to'.
* @return returns TC_CRYPTO_SUCCESS (1)
*         returns TC_CRYPTO_FAIL (0) if:
*                from_len > to_len.
*
* @param to OUT -- destination buffer
* @param to_len IN -- length of destination buffer
* @param from IN -- origin buffer
* @param from_len IN -- length of origin buffer
*/
unsigned int _copy(uint8_t *to, unsigned int to_len,
    const uint8_t *from, unsigned int from_len);

/**
* @brief Set the value 'val' into the buffer 'to', 'len' times.
*
* @param to OUT -- destination buffer
* @param val IN -- value to be set in 'to'
* @param len IN -- number of times the value will be copied
*/
void _set(void *to, uint8_t val, unsigned int len);

/*
* @brief AES specific doubling function, which utilizes
* the finite field used by AES.
* @return Returns a^2
*
* @param a IN/OUT -- value to be doubled
*/
uint8_t _double_byte(uint8_t a);

/*
* @brief Constant-time algorithm to compare if two sequences of bytes are equal
* @return Returns 0 if equal, and non-zero otherwise
*
* @param a IN -- sequence of bytes a
* @param b IN -- sequence of bytes b
* @param size IN -- size of sequences a and b
*/
int _compare(const uint8_t *a, const uint8_t *b, size_t size);

#define TC_SHA256_BLOCK_SIZE (64)
#define TC_SHA256_DIGEST_SIZE (32)
#define TC_SHA256_STATE_BLOCKS (TC_SHA256_DIGEST_SIZE/4)

struct tc_sha256_state_struct {
	unsigned int iv[TC_SHA256_STATE_BLOCKS];
	uint64_t bits_hashed;
	uint8_t leftover[TC_SHA256_BLOCK_SIZE];
	size_t leftover_offset;
};

typedef struct tc_sha256_state_struct *TCSha256State_t;

/**
 *  @brief SHA256 initialization procedure
 *  Initializes s
 *  @return returns TC_CRYPTO_SUCCESS (1)
 *          returns TC_CRYPTO_FAIL (0) if s == NULL
 *  @param s Sha256 state struct
 */
int tc_sha256_init(TCSha256State_t s);

/**
 *  @brief SHA256 update procedure
 *  Hashes data_length bytes addressed by data into state s
 *  @return returns TC_CRYPTO_SUCCESS (1)
 *          returns TC_CRYPTO_FAIL (0) if:
 *                s == NULL,
 *                s->iv == NULL,
 *                data == NULL
 *  @note Assumes s has been initialized by tc_sha256_init
 *  @warning The state buffer 'leftover' is left in memory after processing
 *           If your application intends to have sensitive data in this
 *           buffer, remind to erase it after the data has been processed
 *  @param s Sha256 state struct
 *  @param data message to hash
 *  @param datalen length of message to hash
 */
int tc_sha256_update (TCSha256State_t s, const uint8_t *data, size_t datalen);

/**
 *  @brief SHA256 final procedure
 *  Inserts the completed hash computation into digest
 *  @return returns TC_CRYPTO_SUCCESS (1)
 *          returns TC_CRYPTO_FAIL (0) if:
 *                s == NULL,
 *                s->iv == NULL,
 *                digest == NULL
 *  @note Assumes: s has been initialized by tc_sha256_init
 *        digest points to at least TC_SHA256_DIGEST_SIZE bytes
 *  @warning The state buffer 'leftover' is left in memory after processing
 *           If your application intends to have sensitive data in this
 *           buffer, remind to erase it after the data has been processed
 *  @param digest unsigned eight bit integer
 *  @param Sha256 state struct
 */
int tc_sha256_final(uint8_t *digest, TCSha256State_t s);


void SHA256(uint8_t *digest, const uint8_t *data, size_t datalen, TCSha256State_t s);

#ifdef __cplusplus
}
#endif

#endif /* __TC_SHA256_H__ */
