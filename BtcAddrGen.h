
#ifndef BTCADDRGEN_H
#define BTCADDRGEN_H

#include <inttypes.h>

// --- ECDSA

#define SECP256K1_KEY_SIZE 32

typedef struct Secp256k1Key { 
	unsigned char key[SECP256K1_KEY_SIZE]; 
} Secp256k1Key;

typedef struct Secp256k1UncompressedPublicKey {
	unsigned char x[SECP256K1_KEY_SIZE];
	unsigned char y[SECP256K1_KEY_SIZE];
} Secp256k1UncompressedPublicKey;

int Secp256k1(const Secp256k1Key *privateKey, Secp256k1UncompressedPublicKey *publicKey);

// --- SHA256

#define SHA256_HASH_SIZE 32

typedef struct Sha256Hash {
    unsigned char hash[SHA256_HASH_SIZE];
} Sha256Hash;

void Sha256(const uint8_t *data, uint32_t len, Sha256Hash *hash);

// --- RIPEMD160

#define RIPEMD160_HASH_SIZE 20

typedef struct Ripemd160Hash {
    unsigned char hash[RIPEMD160_HASH_SIZE];
} Ripemd160Hash;

void Ripemd160(const uint8_t *data, uint32_t len, Ripemd160Hash *hash);


// --- BTC

#define BTC_UNCOMPRESSED_PUBLIC_KEY_SIZE ((2* SECP256K1_KEY_SIZE) + 1)
#define BTC_PUBLIC_KEY_SIZE              (    SECP256K1_KEY_SIZE  + 1)
#define BTC_ADDRESS_RAW_SIZE             /   RIPEMD160_HASH_SIZE  + 1)

typedef struct BtcUncompressedPublicKey {
    unsigned char version;
    Secp256k1UncompressedPublicKey key;
} BtcUncompressedPublicKey;

typedef Secp256k1Key BtcPrivateKey;

typedef struct BtcPublicKey {
    unsigned char version;
    Secp256k1Key key;
} BtcPublicKey;

typedef struct BtcAddressRaw {
    unsigned char network;
    Ripemd160Hash hash;
} BtcAddressRaw;

void BtcRaw(const BtcPrivateKey *privateKey, BtcAddressRaw *addressRaw);

#endif // !BTCADDRGEN_H

