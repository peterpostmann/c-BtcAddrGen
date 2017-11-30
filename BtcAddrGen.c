
#include "BtcAddrGen.h"

#include "lib/micro-ecc/uECC.h"
#include "lib/sha256/sha256.h"
#include "lib/ripemd160/rmd160.h"

void vli_print(char *str, uint8_t *vli, unsigned int size) {
    printf("%s ", str);
    for (unsigned i = 0; i<size; ++i) {
        printf("%02X ", (unsigned)vli[i]);
    }
    printf("\n");
}

void BtcRaw(const BtcPrivateKey *privateKey, BtcAddressRaw *addressRaw, int compress)
{
    uint32_t buf[27];

    BtcUncompressedPublicKey pubKey;
    uECC_compute_public_key(privateKey, &(pubKey.key));
    pubKey.version = compress ? (pubKey.key.y[SECP256K1_KEY_SIZE - 1] % 2 == 0 ? 0x02 : 0x03) : 0x04;

    vli_print("Public Key: ", &pubKey, compress ? BTC_PUBLIC_KEY_SIZE : BTC_UNCOMPRESSED_PUBLIC_KEY_SIZE);

    Sha256Hash hash;

    tc_sha256_init(&buf);
    tc_sha256_update(&buf, &pubKey, compress ? BTC_PUBLIC_KEY_SIZE : BTC_UNCOMPRESSED_PUBLIC_KEY_SIZE);
    tc_sha256_final(&hash, &buf);


    vli_print("Sha256: ", &hash, SHA256_HASH_SIZE);

    RMD160Init(&buf);
    RMD160Update(&buf, &hash, SHA256_HASH_SIZE);
    RMD160Final(&(addressRaw->hash), &buf);

    addressRaw->network = 0x00;
}
