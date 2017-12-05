
#include "BtcAddrGen.h"

#include "lib/micro-ecc/uECC.h"
#include "lib/sha256/sha256.h"
#include "lib/ripemd160/rmd160.h"

void vli_print(char *str, uint8_t *vli, unsigned int size) {
    printf("%s ", str);
    for (unsigned i = size; i > 0; i--) {
        printf("%02X ", (unsigned)vli[i - 1]);
    }
    printf("\n");
}

void vli_print_bigEnd(char *str, uint8_t *vli, unsigned int size) {
    printf("%s ", str);
    for (unsigned i = 0; i < size; i++) {
        printf("%02X ", (unsigned)vli[i]);
    }
    printf("\n");
}

void BtcRaw(const BtcPrivateKey *privateKey, BtcAddressRaw *addressRaw)
{
    uint32_t buf[27];

    BtcPublicKey pubKey;
    
    uECC_compute_public_key(privateKey, &pubKey);

    vli_print("Public Key:", &pubKey, BTC_PUBLIC_KEY_SIZE);
       printf("            02 50 86 3A D6 4A 87 AE 8A 2F E8 3C 1A F1 A8 40 3C B5 3F 53 E4 86 D8 51 1D AD 8A 04 88 7E 5B 23 52\n");

    Sha256Hash hash;

    tc_sha256_init(&buf);
    tc_sha256_update(&buf, &pubKey, BTC_PUBLIC_KEY_SIZE);
    tc_sha256_final(&hash, &buf);


    vli_print("Sha256:", &hash, SHA256_HASH_SIZE);
               printf("        0B 7C 28 C9 B7 29 0C 98 D7 43 8E 70 B3 D3 F7 C8 48 FB D7 D1 DC 19 4F F8 3F 4F 7C C9 B1 37 8E 98\n");

    RMD160(&(addressRaw->hash), &hash, SHA256_HASH_SIZE, &buf);

    addressRaw->network = 0x00;
}
