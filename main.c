
#include <stdio.h>

#include "BtcAddrGen.h"

int main()
{
    printf("hello from BtcAddrGen!\n");

	Secp256k1Key privateKey = { .key = { 0x25, 0x17, 0x32, 0x06, 0xA2, 0x29, 0xDB, 0x35, 
                                         0x20, 0x2C, 0x7E, 0xA4, 0xF9, 0xE7, 0x74, 0xE7, 
                                         0xC8, 0xE7, 0x01, 0x47, 0x11, 0xF8, 0x94, 0x6A, 
                                         0x42, 0x7F, 0x30, 0x6A, 0x7B, 0x4A, 0xE1, 0x18} };

    vli_print("Private Key", &privateKey, SECP256K1_KEY_SIZE);
       printf("            18 E1 4A 7B 6A 30 7F 42 6A 94 F8 11 47 01 E7 C8 E7 74 E7 F9 A4 7E 2C 20 35 DB 29 A2 06 32 17 25\n");

    BtcAddressRaw addressRaw;
    BtcRaw(&privateKey, &addressRaw);
    vli_print("Address:", &addressRaw, BTC_ADDRESS_RAW_SIZE);
       printf("         00 F5 4A 58 51 E9 37 2B 87 81 0A 8E 60 CD D2 E7 CF D8 0B 6E 31\n");

    return 0;
}