#if defined(T)

#define D CONCAT3(int, T, _t)

#if defined(STRIPED)
#define NAME parasail_traceback_striped_
#define LOC LOC_STRIPED
#else
#define NAME parasail_traceback_
#define LOC LOC_NOVEC
#endif

static inline void CONCAT(NAME, T) (
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        parasail_result_t *result)
{
    char *q = malloc(sizeof(char)*(lena+lenb));
    char *d = malloc(sizeof(char)*(lena+lenb));
    char *a = malloc(sizeof(char)*(lena+lenb));
    char *qc = q;
    char *dc = d;
    char *ac = a;
    int i = result->end_query;
    int j = result->end_ref;
    int where = PARASAIL_DIAG;
    D *HT = (D*)result->trace_table;
    D *ET = (D*)result->trace_ins_table;
    D *FT = (D*)result->trace_del_table;
#if defined(STRIPED)
    int32_t segWidth = 0;
    int32_t segLen = 0;
    if (result->flag & PARASAIL_FLAG_LANES_1) {
        segWidth = 1;
    }
    if (result->flag & PARASAIL_FLAG_LANES_2) {
        segWidth = 2;
    }
    if (result->flag & PARASAIL_FLAG_LANES_4) {
        segWidth = 4;
    }
    if (result->flag & PARASAIL_FLAG_LANES_8) {
        segWidth = 8;
    }
    if (result->flag & PARASAIL_FLAG_LANES_16) {
        segWidth = 16;
    }
    if (result->flag & PARASAIL_FLAG_LANES_32) {
        segWidth = 32;
    }
    if (result->flag & PARASAIL_FLAG_LANES_64) {
        segWidth = 64;
    }
    segLen = (lena + segWidth - 1) / segWidth;
#endif
    while (i >= 0 && j >= 0) {
        LOC
        assert(i >= 0 && j >= 0);
        if (PARASAIL_DIAG == where) {
            if (HT[loc] == PARASAIL_DIAG) {
                *(qc++) = seqA[i];
                *(dc++) = seqB[j];
                *(ac++) = match_char(seqA[i], seqB[j], matrix);
                --i;
                --j;
            }
            else if (HT[loc] == PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else if (HT[loc] == PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else if (HT[loc] == PARASAIL_ZERO) {
                break;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_INS == where) {
            *(qc++) = '-';
            *(dc++) = seqB[j];
            *(ac++) = ' ';
            --j;
            if (ET[loc] == PARASAIL_DIAG) {
                where = PARASAIL_DIAG;
            }
            else if (ET[loc] == PARASAIL_INS) {
                where = PARASAIL_INS;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_DEL == where) {
            *(qc++) = seqA[i];
            *(dc++) = '-';
            *(ac++) = ' ';
            --i;
            if (FT[loc] == PARASAIL_DIAG) {
                where = PARASAIL_DIAG;
            }
            else if (FT[loc] == PARASAIL_DEL) {
                where = PARASAIL_DEL;
            }
            else {
                assert(0);
            }
        }
        else if (PARASAIL_ZERO == where) {
            break;
        }
        else {
            assert(0);
        }
    }
    *(qc++) = '\0';
    *(dc++) = '\0';
    *(ac++) = '\0';

    if (1) {
        int mch = 0;
        int sim = 0;
        int gap = 0;
        int len = strlen(a);
        printf("Length: %zu\n", strlen(a));
        char *tmp = NULL;
        tmp = parasail_reverse(q, strlen(q));
        printf("%s\n", tmp);
        free(tmp);
        tmp = parasail_reverse(a, strlen(a));
        printf("%s\n", tmp);
        for (i=0; i<len; ++i) {
            if (tmp[i] == '|') { ++mch; ++sim; }
            else if (tmp[i] == ':') ++sim;
            else if (tmp[i] == '.') ;
            else if (tmp[i] == ' ') ++gap;
            else {
                printf("bad char in traceback '%c'\n", tmp[i]);
                assert(0);
            }
        }
        free(tmp);
        tmp = parasail_reverse(d, strlen(d));
        printf("%s\n", tmp);
        free(tmp);
        printf("Length:        %d\n", len);
        printf("Identity:   %d/%d\n", mch, len);
        printf("Similarity: %d/%d\n", sim, len);
        printf("Gaps:       %d/%d\n", gap, len);
    }
    else {
        printf("%s\n", q);
        printf("%s\n", a);
        printf("%s\n", d);
    }

    free(q);
    free(d);
    free(a);
}

#undef D
#undef NAME
#undef LOC

#else

#include "config.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/memory.h"

static inline int weight(
        const char a,
        const char b,
        const parasail_matrix_t *matrix)
{
    return matrix->matrix[matrix->mapper[(unsigned char)a]*matrix->size + matrix->mapper[(unsigned char)b]];
}

static inline char match_char(
        const char a,
        const char b,
        const parasail_matrix_t *matrix)
{
    if (a == b) {
        return '|';
    }
    else {
        int sub = weight(a, b, matrix);
        if (sub > 0) {
            return ':';
        }
        else {
            return '.';
        }
    }

    return 'X'; /* shouldn't happen */
}

#define CONCAT_(X, Y) X##Y
#define CONCAT(X, Y) CONCAT_(X, Y)
#define CONCAT3_(X, Y, Z) X##Y##Z
#define CONCAT3(X, Y, Z) CONCAT3_(X, Y, Z)
#define LOC_NOVEC int loc = i*lenb + j;
#define LOC_STRIPED int loc = j*segLen*segWidth + (i%segLen)*segWidth + (i/segLen);

#define T 8
#include "traceback.c"
#undef T

#define T 8
#define STRIPED
#include "traceback.c"
#undef T
#undef STRIPED

#define T 16
#include "traceback.c"
#undef T

#define T 16
#define STRIPED
#include "traceback.c"
#undef T
#undef STRIPED

#define T 32
#include "traceback.c"
#undef T

#define T 32
#define STRIPED
#include "traceback.c"
#undef T
#undef STRIPED

#define T 64
#include "traceback.c"
#undef T

#define T 64
#define STRIPED
#include "traceback.c"
#undef T
#undef STRIPED

void parasail_traceback(
        const char *seqA,
        int lena,
        const char *seqB,
        int lenb,
        const parasail_matrix_t *matrix,
        parasail_result_t *result)
{
    if (result->flag & PARASAIL_FLAG_STRIPED) {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            parasail_traceback_striped_8(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            parasail_traceback_striped_16(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            parasail_traceback_striped_32(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            parasail_traceback_striped_64(seqA, lena, seqB, lenb, matrix, result);
        }
    }
    else {
        if (result->flag & PARASAIL_FLAG_BITS_8) {
            parasail_traceback_8(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_16) {
            parasail_traceback_16(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_32) {
            parasail_traceback_32(seqA, lena, seqB, lenb, matrix, result);
        }
        else if (result->flag & PARASAIL_FLAG_BITS_64) {
            parasail_traceback_64(seqA, lena, seqB, lenb, matrix, result);
        }
    }
}

#endif

