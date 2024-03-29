/* RMD160.H - header file for RMD160.C
 */
#ifndef _RMD160_H_
#define _RMD160_H_

#include <sys/types.h>

#define RMD160_BLOCKBYTES 64
#define RMD160_BLOCKWORDS 16

#define RMD160_HASHBYTES 20
#define RMD160_HASHWORDS 5

/* For compatibility */
#define RIPEMD160_BLOCKBYTES 64
#define RIPEMD160_BLOCKWORDS 16

#define RIPEMD160_HASHBYTES 20
#define RIPEMD160_HASHWORDS 5

/* RIPEMD160 context. */
typedef struct RMD160Context {
 u_int32_t key[RIPEMD160_BLOCKWORDS];
 u_int32_t iv[RIPEMD160_HASHWORDS];
 u_int32_t bytesHi, bytesLo;
} RMD160_CTX;

#define RIPEMD160Context RMD160_CTX

#include <sys/cdefs.h>

__BEGIN_DECLS
void RMD160(unsigned char[RMD160_HASHBYTES], const u_int8_t *, size_t, RMD160_CTX *);
__END_DECLS

#endif /* _RMD160_H_ */
