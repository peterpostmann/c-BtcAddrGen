
#include "BtcAddrGen.h"


// --- ECDSA

#include "lib/micro-ecc/uECC.h"

#if SECP256K1_KEY_SIZE != uECC_BYTES
#error "Configuration Error: SECP256K1_KEY_SIZE must match uECC_BYTES"
#endif

int Secp256k1(const Secp256k1Key *privateKey, Secp256k1UncompressedPublicKey *publicKey)
{	
	return uECC_compute_public_key(privateKey, publicKey);
}


// --- SHA256

#include "lib/sha256/sha256.h"

void Sha256(const uint8_t *data, uint32_t len, Sha256Hash *hash)
{
    sha256_context ctx;

    sha256_init(&ctx);
    sha256_hash(&ctx, data, len);
    sha256_done(&ctx, hash);
}


// --- RIPEMD160

#include "lib/ripemd160/rmd160.h"

void Ripemd160(const uint8_t *data, uint32_t len, Ripemd160Hash *hash)
{
    RIPEMD160Context ctx;

    RMD160Init(&ctx);
    RMD160Update(&ctx, data, len);
    RMD160Final(hash, &ctx);
}
