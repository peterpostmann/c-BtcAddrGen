
#include "BtcAddrGen.h"

#include "lib/micro-ecc/uECC.h"
#include "lib/sha256/sha256.h"
#include "lib/ripemd160/rmd160.h"

void BtcRaw(const BtcPrivateKey *privateKey, BtcAddressRaw *addressRaw)
{
    uint32_t buf[26];

    BtcUncompressedPublicKey pubKey = { .version = 0x04 };
    uECC_compute_public_key(privateKey, &(pubKey.key));

    Sha256Hash hash;

    sha256_init(&buf);
    sha256_hash(&buf, &pubKey, BTC_UNCOMPRESSED_PUBLIC_KEY_SIZE);
    sha256_done(&buf, &hash);

    RMD160Init(&buf);
    RMD160Update(&buf, &hash, SHA256_HASH_SIZE);
    RMD160Final(addressRaw, &buf);

    addressRaw->network = 0x00;
}
