/* Copyright 2014, Kenneth MacKay. Licensed under the BSD 2-clause license. */

#include "uECC.h"

#if __STDC_VERSION__ >= 199901L
    #define RESTRICT restrict
#else
    #define RESTRICT
#endif

#define uECC_BYTES 32

#define MAX_TRIES 64

typedef uint64_t uECC_word_t;
typedef unsigned __int128 uECC_dword_t;
typedef unsigned wordcount_t;
typedef int swordcount_t;
typedef int bitcount_t;
typedef int cmpresult_t;

#define HIGH_BIT_SET 0x8000000000000000ull
#define uECC_WORD_BITS 64
#define uECC_WORD_BITS_SHIFT 6
#define uECC_WORD_BITS_MASK 0x03F

#define Curve_P {0xFFFFFFFEFFFFFC2Full, 0xFFFFFFFFFFFFFFFFull, \
                 0xFFFFFFFFFFFFFFFFull, 0xFFFFFFFFFFFFFFFFull}

#define Curve_B {0x0000000000000007ull, 0x0000000000000000ull, \
                 0x0000000000000000ull, 0x0000000000000000ull}

#define Curve_G { \
    {0x59F2815B16F81798ull, 0x029BFCDB2DCE28D9ull, 0x55A06295CE870B07ull, 0x79BE667EF9DCBBACull}, \
    {0x9C47D08FFB10D4B8ull, 0xFD17B448A6855419ull, 0x5DA4FBFC0E1108A8ull, 0x483ADA7726A3C465ull}}

#define Curve_N {0xBFD25E8CD0364141ull, 0xBAAEDCE6AF48A03Bull, \
                 0xFFFFFFFFFFFFFFFEull, 0xFFFFFFFFFFFFFFFFull}

#define uECC_WORDS 4
#define uECC_N_WORDS 4

typedef struct EccPoint {
    uECC_word_t x[uECC_WORDS];
    uECC_word_t y[uECC_WORDS];
} EccPoint;

static const uECC_word_t curve_p[uECC_WORDS] = Curve_P;
static const uECC_word_t curve_b[uECC_WORDS] = Curve_B;
static const EccPoint curve_G = Curve_G;
static const uECC_word_t curve_n[uECC_N_WORDS] = Curve_N;
static void vli_clear(uECC_word_t *vli);
static uECC_word_t vli_isZero(const uECC_word_t *vli);
static uECC_word_t vli_testBit(const uECC_word_t *vli, bitcount_t bit);
static bitcount_t vli_numBits(const uECC_word_t *vli, wordcount_t max_words);
static void vli_set(uECC_word_t *dest, const uECC_word_t *src);
static cmpresult_t vli_cmp(const uECC_word_t *left, const uECC_word_t *right);
static cmpresult_t vli_equal(const uECC_word_t *left, const uECC_word_t *right);
static void vli_rshift1(uECC_word_t *vli);
static uECC_word_t vli_add(uECC_word_t *result,
                           const uECC_word_t *left,
                           const uECC_word_t *right);
static uECC_word_t vli_sub(uECC_word_t *result,
                           const uECC_word_t *left,
                           const uECC_word_t *right);
static void vli_mult(uECC_word_t *result, const uECC_word_t *left, const uECC_word_t *right);
static void vli_modAdd(uECC_word_t *result,
                       const uECC_word_t *left,
                       const uECC_word_t *right,
                       const uECC_word_t *mod);
static void vli_modSub(uECC_word_t *result,
                       const uECC_word_t *left,
                       const uECC_word_t *right,
                       const uECC_word_t *mod);
static void vli_mmod_fast(uECC_word_t *RESTRICT result, uECC_word_t *RESTRICT product);
static void vli_modMult_fast(uECC_word_t *result,
                             const uECC_word_t *left,
                             const uECC_word_t *right);
static void vli_modInv(uECC_word_t *result, const uECC_word_t *input, const uECC_word_t *mod);
#if uECC_SQUARE_FUNC
static void vli_square(uECC_word_t *result, const uECC_word_t *left);
static void vli_modSquare_fast(uECC_word_t *result, const uECC_word_t *left);
#endif

static void vli_clear(uECC_word_t *vli) {
    wordcount_t i;
    for (i = 0; i < uECC_WORDS; ++i) {
        vli[i] = 0;
    }
}

/* Returns 1 if vli == 0, 0 otherwise. */
static uECC_word_t vli_isZero(const uECC_word_t *vli) {
    wordcount_t i;
    for (i = 0; i < uECC_WORDS; ++i) {
        if (vli[i]) {
            return 0;
        }
    }
    return 1;
}

/* Returns nonzero if bit 'bit' of vli is set. */
static uECC_word_t vli_testBit(const uECC_word_t *vli, bitcount_t bit) {
    return (vli[bit >> uECC_WORD_BITS_SHIFT] & ((uECC_word_t)1 << (bit & uECC_WORD_BITS_MASK)));
}

/* Counts the number of words in vli. */
static wordcount_t vli_numDigits(const uECC_word_t *vli, wordcount_t max_words) {
    swordcount_t i;
    /* Search from the end until we find a non-zero digit.
       We do it in reverse because we expect that most digits will be nonzero. */
    for (i = max_words - 1; i >= 0 && vli[i] == 0; --i) {
    }

    return (i + 1);
}

/* Counts the number of bits required to represent vli. */
static bitcount_t vli_numBits(const uECC_word_t *vli, wordcount_t max_words) {
    uECC_word_t i;
    uECC_word_t digit;

    wordcount_t num_digits = vli_numDigits(vli, max_words);
    if (num_digits == 0) {
        return 0;
    }

    digit = vli[num_digits - 1];
    for (i = 0; digit; ++i) {
        digit >>= 1;
    }

    return (((bitcount_t)(num_digits - 1) << uECC_WORD_BITS_SHIFT) + i);
}

/* Sets dest = src. */
static void vli_set(uECC_word_t *dest, const uECC_word_t *src) {
    wordcount_t i;
    for (i = 0; i < uECC_WORDS; ++i) {
        dest[i] = src[i];
    }
}

/* Returns sign of left - right. */
static cmpresult_t vli_cmp(const uECC_word_t *left, const uECC_word_t *right) {
    swordcount_t i;
    for (i = uECC_WORDS - 1; i >= 0; --i) {
        if (left[i] > right[i]) {
            return 1;
        } else if (left[i] < right[i]) {
            return -1;
        }
    }
    return 0;
}

static cmpresult_t vli_equal(const uECC_word_t *left, const uECC_word_t *right) {
    uECC_word_t result = 0;
    swordcount_t i;
    for (i = uECC_WORDS - 1; i >= 0; --i) {
        result |= (left[i] ^ right[i]);
    }
    return (result == 0);
}

/* Computes vli = vli >> 1. */
static void vli_rshift1(uECC_word_t *vli) {
    uECC_word_t *end = vli;
    uECC_word_t carry = 0;

    vli += uECC_WORDS;
    while (vli-- > end) {
        uECC_word_t temp = *vli;
        *vli = (temp >> 1) | carry;
        carry = temp << (uECC_WORD_BITS - 1);
    }
}

/* Computes result = left + right, returning carry. Can modify in place. */
static uECC_word_t vli_add(uECC_word_t *result, const uECC_word_t *left, const uECC_word_t *right) {
    uECC_word_t carry = 0;
    wordcount_t i;
    for (i = 0; i < uECC_WORDS; ++i) {
        uECC_word_t sum = left[i] + right[i] + carry;
        if (sum != left[i]) {
            carry = (sum < left[i]);
        }
        result[i] = sum;
    }
    return carry;
}

/* Computes result = left - right, returning borrow. Can modify in place. */
static uECC_word_t vli_sub(uECC_word_t *result, const uECC_word_t *left, const uECC_word_t *right) {
    uECC_word_t borrow = 0;
    wordcount_t i;
    for (i = 0; i < uECC_WORDS; ++i) {
        uECC_word_t diff = left[i] - right[i] - borrow;
        if (diff != left[i]) {
            borrow = (diff > left[i]);
        }
        result[i] = diff;
    }
    return borrow;
}

static void muladd(uECC_word_t a,
                   uECC_word_t b,
                   uECC_word_t *r0,
                   uECC_word_t *r1,
                   uECC_word_t *r2) {
    uECC_dword_t p = (uECC_dword_t)a * b;
    uECC_dword_t r01 = ((uECC_dword_t)(*r1) << uECC_WORD_BITS) | *r0;
    r01 += p;
    *r2 += (r01 < p);
    *r1 = r01 >> uECC_WORD_BITS;
    *r0 = (uECC_word_t)r01;

}

static void vli_mult(uECC_word_t *result, const uECC_word_t *left, const uECC_word_t *right) {
    uECC_word_t r0 = 0;
    uECC_word_t r1 = 0;
    uECC_word_t r2 = 0;
    wordcount_t i, k;

    /* Compute each digit of result in sequence, maintaining the carries. */
    for (k = 0; k < uECC_WORDS; ++k) {
        for (i = 0; i <= k; ++i) {
            muladd(left[i], right[k - i], &r0, &r1, &r2);
        }
        result[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }
    for (k = uECC_WORDS; k < uECC_WORDS * 2 - 1; ++k) {
        for (i = (k + 1) - uECC_WORDS; i < uECC_WORDS; ++i) {
            muladd(left[i], right[k - i], &r0, &r1, &r2);
        }
        result[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }
    result[uECC_WORDS * 2 - 1] = r0;
}

#if uECC_SQUARE_FUNC

static void mul2add(uECC_word_t a,
                    uECC_word_t b,
                    uECC_word_t *r0,
                    uECC_word_t *r1,
                    uECC_word_t *r2) {
    uECC_dword_t p = (uECC_dword_t)a * b;
    uECC_dword_t r01 = ((uECC_dword_t)(*r1) << uECC_WORD_BITS) | *r0;
    *r2 += (p >> (uECC_WORD_BITS * 2 - 1));
    p *= 2;
    r01 += p;
    *r2 += (r01 < p);
    *r1 = r01 >> uECC_WORD_BITS;
    *r0 = (uECC_word_t)r01;
}

static void vli_square(uECC_word_t *result, const uECC_word_t *left) {
    uECC_word_t r0 = 0;
    uECC_word_t r1 = 0;
    uECC_word_t r2 = 0;

    wordcount_t i, k;

    for (k = 0; k < uECC_WORDS * 2 - 1; ++k) {
        uECC_word_t min = (k < uECC_WORDS ? 0 : (k + 1) - uECC_WORDS);
        for (i = min; i <= k && i <= k - i; ++i) {
            if (i < k-i) {
                mul2add(left[i], left[k - i], &r0, &r1, &r2);
            } else {
                muladd(left[i], left[k - i], &r0, &r1, &r2);
            }
        }
        result[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }

    result[uECC_WORDS * 2 - 1] = r0;
}

#else /* uECC_SQUARE_FUNC */

#define vli_square(result, left, size) vli_mult((result), (left), (left), (size))

#endif /* uECC_SQUARE_FUNC */


/* Computes result = (left + right) % mod.
   Assumes that left < mod and right < mod, and that result does not overlap mod. */
static void vli_modAdd(uECC_word_t *result,
                       const uECC_word_t *left,
                       const uECC_word_t *right,
                       const uECC_word_t *mod) {
    uECC_word_t carry = vli_add(result, left, right);
    if (carry || vli_cmp(result, mod) >= 0) {
        /* result > mod (result = mod + remainder), so subtract mod to get remainder. */
        vli_sub(result, result, mod);
    }
}

/* Computes result = (left - right) % mod.
   Assumes that left < mod and right < mod, and that result does not overlap mod. */
static void vli_modSub(uECC_word_t *result,
                       const uECC_word_t *left,
                       const uECC_word_t *right,
                       const uECC_word_t *mod) {
    uECC_word_t l_borrow = vli_sub(result, left, right);
    if (l_borrow) {
        /* In this case, result == -diff == (max int) - diff. Since -x % d == d - x,
           we can get the correct result from result + mod (with overflow). */
        vli_add(result, result, mod);
    }
}

#define vli_modSub_fast(result, left, right) vli_modSub((result), (left), (right), curve_p)

static void omega_mult(uECC_word_t * RESTRICT result, const uECC_word_t * RESTRICT right);

/* Computes result = product % curve_p
    see http://www.isys.uni-klu.ac.at/PDF/2001-0126-MT.pdf page 354

    Note that this only works if log2(omega) < log2(p) / 2 */
static void vli_mmod_fast(uECC_word_t *RESTRICT result, uECC_word_t *RESTRICT product) {
    uECC_word_t tmp[2 * uECC_WORDS];
    uECC_word_t carry;

    vli_clear(tmp);
    vli_clear(tmp + uECC_WORDS);

    omega_mult(tmp, product + uECC_WORDS); /* (Rq, q) = q * c */

    carry = vli_add(result, product, tmp); /* (C, r) = r + q       */
    vli_clear(product);
    omega_mult(product, tmp + uECC_WORDS); /* Rq*c */
    carry += vli_add(result, result, product); /* (C1, r) = r + Rq*c */

    while (carry > 0) {
        --carry;
        vli_sub(result, result, curve_p);
    }
    if (vli_cmp(result, curve_p) > 0) {
        vli_sub(result, result, curve_p);
    }
}

static void omega_mult(uint64_t * RESTRICT result, const uint64_t * RESTRICT right) {
    uECC_word_t r0 = 0;
    uECC_word_t r1 = 0;
    uECC_word_t r2 = 0;
    wordcount_t k;

    /* Multiply by (2^32 + 2^9 + 2^8 + 2^7 + 2^6 + 2^4 + 1). */
    for (k = 0; k < uECC_WORDS; ++k) {
        muladd(0x1000003D1ull, right[k], &r0, &r1, &r2);
        result[k] = r0;
        r0 = r1;
        r1 = r2;
        r2 = 0;
    }
    result[uECC_WORDS] = r0;
}


/* Computes result = (left * right) % curve_p. */
static void vli_modMult_fast(uECC_word_t *result,
                             const uECC_word_t *left,
                             const uECC_word_t *right) {
    uECC_word_t product[2 * uECC_WORDS];
    vli_mult(product, left, right);
    vli_mmod_fast(result, product);
}

#if uECC_SQUARE_FUNC

/* Computes result = left^2 % curve_p. */
static void vli_modSquare_fast(uECC_word_t *result, const uECC_word_t *left) {
    uECC_word_t product[2 * uECC_WORDS];
    vli_square(product, left);
    vli_mmod_fast(result, product);
}

#else /* uECC_SQUARE_FUNC */

#define vli_modSquare_fast(result, left) vli_modMult_fast((result), (left), (left))

#endif /* uECC_SQUARE_FUNC */


#define EVEN(vli) (!(vli[0] & 1))
/* Computes result = (1 / input) % mod. All VLIs are the same size.
   See "From Euclid's GCD to Montgomery Multiplication to the Great Divide"
   https://labs.oracle.com/techrep/2001/smli_tr-2001-95.pdf */
static void vli_modInv(uECC_word_t *result, const uECC_word_t *input, const uECC_word_t *mod) {
    uECC_word_t a[uECC_WORDS], b[uECC_WORDS], u[uECC_WORDS], v[uECC_WORDS];
    uECC_word_t carry;
    cmpresult_t cmpResult;

    if (vli_isZero(input)) {
        vli_clear(result);
        return;
    }

    vli_set(a, input);
    vli_set(b, mod);
    vli_clear(u);
    u[0] = 1;
    vli_clear(v);
    while ((cmpResult = vli_cmp(a, b)) != 0) {
        carry = 0;
        if (EVEN(a)) {
            vli_rshift1(a);
            if (!EVEN(u)) {
                carry = vli_add(u, u, mod);
            }
            vli_rshift1(u);
            if (carry) {
                u[uECC_WORDS - 1] |= HIGH_BIT_SET;
            }
        } else if (EVEN(b)) {
            vli_rshift1(b);
            if (!EVEN(v)) {
                carry = vli_add(v, v, mod);
            }
            vli_rshift1(v);
            if (carry) {
                v[uECC_WORDS - 1] |= HIGH_BIT_SET;
            }
        } else if (cmpResult > 0) {
            vli_sub(a, a, b);
            vli_rshift1(a);
            if (vli_cmp(u, v) < 0) {
                vli_add(u, u, mod);
            }
            vli_sub(u, u, v);
            if (!EVEN(u)) {
                carry = vli_add(u, u, mod);
            }
            vli_rshift1(u);
            if (carry) {
                u[uECC_WORDS - 1] |= HIGH_BIT_SET;
            }
        } else {
            vli_sub(b, b, a);
            vli_rshift1(b);
            if (vli_cmp(v, u) < 0) {
                vli_add(v, v, mod);
            }
            vli_sub(v, v, u);
            if (!EVEN(v)) {
                carry = vli_add(v, v, mod);
            }
            vli_rshift1(v);
            if (carry) {
                v[uECC_WORDS - 1] |= HIGH_BIT_SET;
            }
        }
    }
    vli_set(result, u);
}

/* ------ Point operations ------ */

/* Point multiplication algorithm using Montgomery's ladder with co-Z coordinates.
From http://eprint.iacr.org/2011/338.pdf
*/

/* Double in place */
static void EccPoint_double_jacobian(uECC_word_t * RESTRICT X1,
                                     uECC_word_t * RESTRICT Y1,
                                     uECC_word_t * RESTRICT Z1) {
    /* t1 = X, t2 = Y, t3 = Z */
    uECC_word_t t4[uECC_WORDS];
    uECC_word_t t5[uECC_WORDS];

    if (vli_isZero(Z1)) {
        return;
    }

    vli_modSquare_fast(t5, Y1);   /* t5 = y1^2 */
    vli_modMult_fast(t4, X1, t5); /* t4 = x1*y1^2 = A */
    vli_modSquare_fast(X1, X1);   /* t1 = x1^2 */
    vli_modSquare_fast(t5, t5);   /* t5 = y1^4 */
    vli_modMult_fast(Z1, Y1, Z1); /* t3 = y1*z1 = z3 */

    vli_modAdd(Y1, X1, X1, curve_p); /* t2 = 2*x1^2 */
    vli_modAdd(Y1, Y1, X1, curve_p); /* t2 = 3*x1^2 */
    if (vli_testBit(Y1, 0)) {
        uECC_word_t carry = vli_add(Y1, Y1, curve_p);
        vli_rshift1(Y1);
        Y1[uECC_WORDS - 1] |= carry << (uECC_WORD_BITS - 1);
    } else {
        vli_rshift1(Y1);
    }
    /* t2 = 3/2*(x1^2) = B */

    vli_modSquare_fast(X1, Y1);      /* t1 = B^2 */
    vli_modSub(X1, X1, t4, curve_p); /* t1 = B^2 - A */
    vli_modSub(X1, X1, t4, curve_p); /* t1 = B^2 - 2A = x3 */

    vli_modSub(t4, t4, X1, curve_p); /* t4 = A - x3 */
    vli_modMult_fast(Y1, Y1, t4);    /* t2 = B * (A - x3) */
    vli_modSub(Y1, Y1, t5, curve_p); /* t2 = B * (A - x3) - y1^4 = y3 */
}

/* Modify (x1, y1) => (x1 * z^2, y1 * z^3) */
static void apply_z(uECC_word_t * RESTRICT X1,
                    uECC_word_t * RESTRICT Y1,
                    const uECC_word_t * RESTRICT Z) {
    uECC_word_t t1[uECC_WORDS];

    vli_modSquare_fast(t1, Z);    /* z^2 */
    vli_modMult_fast(X1, X1, t1); /* x1 * z^2 */
    vli_modMult_fast(t1, t1, Z);  /* z^3 */
    vli_modMult_fast(Y1, Y1, t1); /* y1 * z^3 */
}

/* P = (x1, y1) => 2P, (x2, y2) => P' */
static void XYcZ_initial_double(uECC_word_t * RESTRICT X1,
                                uECC_word_t * RESTRICT Y1,
                                uECC_word_t * RESTRICT X2,
                                uECC_word_t * RESTRICT Y2,
                                const uECC_word_t * RESTRICT initial_Z) {
    uECC_word_t z[uECC_WORDS];
    if (initial_Z) {
        vli_set(z, initial_Z);
    } else {
        vli_clear(z);
        z[0] = 1;
    }

    vli_set(X2, X1);
    vli_set(Y2, Y1);

    apply_z(X1, Y1, z);
    EccPoint_double_jacobian(X1, Y1, z);
    apply_z(X2, Y2, z);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
   Output P' = (x1', y1', Z3), P + Q = (x3, y3, Z3)
   or P => P', Q => P + Q
*/
static void XYcZ_add(uECC_word_t * RESTRICT X1,
                     uECC_word_t * RESTRICT Y1,
                     uECC_word_t * RESTRICT X2,
                     uECC_word_t * RESTRICT Y2) {
    /* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
    uECC_word_t t5[uECC_WORDS];

    vli_modSub_fast(t5, X2, X1);  /* t5 = x2 - x1 */
    vli_modSquare_fast(t5, t5);   /* t5 = (x2 - x1)^2 = A */
    vli_modMult_fast(X1, X1, t5); /* t1 = x1*A = B */
    vli_modMult_fast(X2, X2, t5); /* t3 = x2*A = C */
    vli_modSub_fast(Y2, Y2, Y1);  /* t4 = y2 - y1 */
    vli_modSquare_fast(t5, Y2);   /* t5 = (y2 - y1)^2 = D */

    vli_modSub_fast(t5, t5, X1);  /* t5 = D - B */
    vli_modSub_fast(t5, t5, X2);  /* t5 = D - B - C = x3 */
    vli_modSub_fast(X2, X2, X1);  /* t3 = C - B */
    vli_modMult_fast(Y1, Y1, X2); /* t2 = y1*(C - B) */
    vli_modSub_fast(X2, X1, t5);  /* t3 = B - x3 */
    vli_modMult_fast(Y2, Y2, X2); /* t4 = (y2 - y1)*(B - x3) */
    vli_modSub_fast(Y2, Y2, Y1);  /* t4 = y3 */

    vli_set(X2, t5);
}

/* Input P = (x1, y1, Z), Q = (x2, y2, Z)
   Output P + Q = (x3, y3, Z3), P - Q = (x3', y3', Z3)
   or P => P - Q, Q => P + Q
*/
static void XYcZ_addC(uECC_word_t * RESTRICT X1,
                      uECC_word_t * RESTRICT Y1,
                      uECC_word_t * RESTRICT X2,
                      uECC_word_t * RESTRICT Y2) {
    /* t1 = X1, t2 = Y1, t3 = X2, t4 = Y2 */
    uECC_word_t t5[uECC_WORDS];
    uECC_word_t t6[uECC_WORDS];
    uECC_word_t t7[uECC_WORDS];

    vli_modSub_fast(t5, X2, X1);     /* t5 = x2 - x1 */
    vli_modSquare_fast(t5, t5);      /* t5 = (x2 - x1)^2 = A */
    vli_modMult_fast(X1, X1, t5);    /* t1 = x1*A = B */
    vli_modMult_fast(X2, X2, t5);    /* t3 = x2*A = C */
    vli_modAdd(t5, Y2, Y1, curve_p); /* t5 = y2 + y1 */
    vli_modSub_fast(Y2, Y2, Y1);     /* t4 = y2 - y1 */

    vli_modSub_fast(t6, X2, X1);     /* t6 = C - B */
    vli_modMult_fast(Y1, Y1, t6);    /* t2 = y1 * (C - B) = E */
    vli_modAdd(t6, X1, X2, curve_p); /* t6 = B + C */
    vli_modSquare_fast(X2, Y2);      /* t3 = (y2 - y1)^2 = D */
    vli_modSub_fast(X2, X2, t6);     /* t3 = D - (B + C) = x3 */

    vli_modSub_fast(t7, X1, X2);  /* t7 = B - x3 */
    vli_modMult_fast(Y2, Y2, t7); /* t4 = (y2 - y1)*(B - x3) */
    vli_modSub_fast(Y2, Y2, Y1);  /* t4 = (y2 - y1)*(B - x3) - E = y3 */

    vli_modSquare_fast(t7, t5);   /* t7 = (y2 + y1)^2 = F */
    vli_modSub_fast(t7, t7, t6);  /* t7 = F - (B + C) = x3' */
    vli_modSub_fast(t6, t7, X1);  /* t6 = x3' - B */
    vli_modMult_fast(t6, t6, t5); /* t6 = (y2 + y1)*(x3' - B) */
    vli_modSub_fast(Y1, t6, Y1);  /* t2 = (y2 + y1)*(x3' - B) - E = y3' */

    vli_set(X1, t7);
}

int uECC_compute_public_key(const uint8_t private_key[uECC_BYTES],
    uint8_t public_key[uECC_BYTES * 2]) {
    uECC_word_t private[uECC_WORDS];
    EccPoint public;

    EccPoint * RESTRICT result = &public;
    EccPoint * RESTRICT point = &curve_G;
    uECC_word_t * RESTRICT scalar = private_key;
    uECC_word_t * RESTRICT initialZ = 0;
    bitcount_t numBits = vli_numBits(private_key, uECC_WORDS);

    uECC_word_t Rx[2][uECC_WORDS];
    uECC_word_t Ry[2][uECC_WORDS];
    uECC_word_t z[uECC_WORDS];
    bitcount_t i;
    uECC_word_t nb;

    vli_set(Rx[1], point->x);
    vli_set(Ry[1], point->y);

    XYcZ_initial_double(Rx[1], Ry[1], Rx[0], Ry[0], initialZ);

    for (i = numBits - 2; i > 0; --i) {
        nb = !vli_testBit(scalar, i);
        XYcZ_addC(Rx[1 - nb], Ry[1 - nb], Rx[nb], Ry[nb]);
        XYcZ_add(Rx[nb], Ry[nb], Rx[1 - nb], Ry[1 - nb]);
    }

    nb = !vli_testBit(scalar, 0);
    XYcZ_addC(Rx[1 - nb], Ry[1 - nb], Rx[nb], Ry[nb]);

    /* Find final 1/Z value. */
    vli_modSub_fast(z, Rx[1], Rx[0]);   /* X1 - X0 */
    vli_modMult_fast(z, z, Ry[1 - nb]); /* Yb * (X1 - X0) */
    vli_modMult_fast(z, z, point->x); /* xP * Yb * (X1 - X0) */
    vli_modInv(z, z, curve_p);          /* 1 / (xP * Yb * (X1 - X0)) */
    vli_modMult_fast(z, z, point->y); /* yP / (xP * Yb * (X1 - X0)) */
    vli_modMult_fast(z, z, Rx[1 - nb]); /* Xb * yP / (xP * Yb * (X1 - X0)) */
                                        /* End 1/Z calculation */

    XYcZ_add(Rx[nb], Ry[nb], Rx[1 - nb], Ry[1 - nb]);
    apply_z(Rx[0], Ry[0], z);

    vli_set(public_key, Rx[0]);
    *(public_key + uECC_BYTES) = Ry[0][uECC_WORDS - 1] % 2 ? 0x02 : 0x03;
}
