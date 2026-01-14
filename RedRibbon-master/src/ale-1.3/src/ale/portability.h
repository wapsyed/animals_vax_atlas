#ifndef __PORTABILITY_H
#define __PORTABILITY_H

#include <limits.h>
// #include <netdb.h>        // REMOVIDO: Nao existe no Windows
#include <sys/types.h>    // MANTIDO: Necessario
#include <stdint.h>       // MANTIDO: Necessario
#include <math.h>         // MANTIDO: Necessario
// #include <sys/random.h>   // REMOVIDO: Linux only

//  // COMENTADO: Erro de script

#ifndef MAX
# define MAX(a,b) ( ((a) < (b)) ? (b) : (a))
#endif

#ifndef MIN
# define MIN(a,b) ( ((a) <= (b)) ? (a) : (b))
#endif

#ifndef PATH_MAX
# define PATH_MAX 4096
#endif

#ifndef NI_MAXHOST
#define NI_MAXHOST 1025
#endif

/* --- INICIO DA CORRECAO PARA WINDOWS --- */
#if defined(_WIN32) || defined(_WIN64)
   // No Windows (MinGW/Rtools), usamos __builtin_bswap64 para inverter bytes
   // se estivermos em little-endian (padrao Intel/AMD)
   #define htobe64(x) __builtin_bswap64(x)
   #define be64toh(x) __builtin_bswap64(x)
   
#elif defined(__linux__)
#  include <endian.h>
#elif defined(__FreeBSD__) || defined(__NetBSD__)
#  include <sys/endian.h>
#elif defined(__OpenBSD__)
#  include <sys/types.h>
#  define be16toh(x) betoh16(x)
#  define be32toh(x) betoh32(x)
#  define be64toh(x) betoh64(x)
#endif
/* --- FIM DA CORRECAO --- */

#ifndef htonll
#define htonll(x) htobe64(x)
#endif

#ifndef ntohll
#define ntohll(x) be64toh(x)
#endif

#ifndef LOG_PERROR
# define LOG_PERROR 0
#endif

#ifndef HAVE_POSIX_FADVISE
#  define posix_fadvise portability_posix_fadvise
#  ifndef POSIX_FADV_SEQUENTIAL
#    define POSIX_FADV_SEQUENTIAL
