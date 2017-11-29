
#include <stdio.h>

#include "BtcAddrGen.h"

void vli_print(char *str, uint8_t *vli, unsigned int size) {
	printf("%s ", str);
	for (unsigned i = 0; i<size; ++i) {
		printf("%02X ", (unsigned)vli[i]);
	}
	printf("\n");
}

int main()
{
    printf("hello from BtcAddrGen!\n");

	Secp256k1Key privateKey = { .key = { 0x18, 0xE1, 0x4A, 0x7B, 0x6A, 0x30, 0x7F, 0x42, 0x6A, 0x94, 0xF8, 0x11, 0x47, 0x01, 0xE7, 0xC8, 0xE7, 0x74, 0xE7, 0xF9, 0xA4, 0x7E, 0x2C, 0x20, 0x35, 0xDB, 0x29, 0xA2, 0x06, 0x32, 0x17, 0x25 } };
	vli_print("Private Key: ", &privateKey, SECP256K1_KEY_SIZE);
    
    BtcAddressRaw rawAddress;
    BtcRaw(&privateKey, &rawAddress);
    vli_print("Address: ", &rawAddress, RIPEMD160_HASH_SIZE);

    return 0;
}