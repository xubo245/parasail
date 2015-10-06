/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdint.h>
#include <stdlib.h>

#include <emmintrin.h>
#include <smmintrin.h>

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/internal_sse.h"

#define SWAP(A,B) { __m128i* tmp = A; A = B; B = tmp; }
#define SWAP3(A,B,C) { __m128i* tmp = A; A = B; B = C; C = tmp; }

#define NEG_INF (INT64_MIN/(int64_t)(2))

static inline __m128i _mm_cmpgt_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? 0xFFFFFFFFFFFFFFFF : 0;
    A.v[1] = (A.v[1]>B.v[1]) ? 0xFFFFFFFFFFFFFFFF : 0;
    return A.m;
}

static inline __m128i _mm_max_epi64_rpl(__m128i a, __m128i b) {
    __m128i_64_t A;
    __m128i_64_t B;
    A.m = a;
    B.m = b;
    A.v[0] = (A.v[0]>B.v[0]) ? A.v[0] : B.v[0];
    A.v[1] = (A.v[1]>B.v[1]) ? A.v[1] : B.v[1];
    return A.m;
}

static inline int64_t _mm_hmax_epi64_rpl(__m128i a) {
    a = _mm_max_epi64_rpl(a, _mm_srli_si128(a, 8));
    return _mm_extract_epi64(a, 0);
}


#if defined(PARASAIL_TABLE) || defined(PARASAIL_TRACE)
static inline void arr_store(
        int *array,
        __m128i vH,
        int32_t t,
        int32_t seglen,
        int32_t d,
        int32_t dlen)
{
    array[(0*seglen+t)*dlen + d] = (int64_t)_mm_extract_epi64(vH, 0);
    array[(1*seglen+t)*dlen + d] = (int64_t)_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_ROWCOL
static inline void arr_store_col(
        int *col,
        __m128i vH,
        int32_t t,
        int32_t seglen)
{
    col[0*seglen+t] = (int64_t)_mm_extract_epi64(vH, 0);
    col[1*seglen+t] = (int64_t)_mm_extract_epi64(vH, 1);
}
#endif

#ifdef PARASAIL_TABLE
#define FNAME parasail_sw_table_striped_sse41_128_64
#define PNAME parasail_sw_table_striped_profile_sse41_128_64
#else
#ifdef PARASAIL_ROWCOL
#define FNAME parasail_sw_rowcol_striped_sse41_128_64
#define PNAME parasail_sw_rowcol_striped_profile_sse41_128_64
#else
#ifdef PARASAIL_TRACE
#define FNAME parasail_sw_trace_striped_sse41_128_64
#define PNAME parasail_sw_trace_striped_profile_sse41_128_64
#else
#define FNAME parasail_sw_striped_sse41_128_64
#define PNAME parasail_sw_striped_profile_sse41_128_64
#endif
#endif
#endif

parasail_result_t* FNAME(
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_profile_t *profile = parasail_profile_create_sse_128_64(s1, s1Len, matrix);
    parasail_result_t *result = PNAME(profile, s2, s2Len, open, gap);
    parasail_profile_free(profile);
    return result;
}

parasail_result_t* PNAME(
        const parasail_profile_t * const restrict profile,
        const char * const restrict s2, const int s2Len,
        const int open, const int gap)
{
    int32_t i = 0;
    int32_t j = 0;
    int32_t k = 0;
    int32_t end_query = 0;
    int32_t end_ref = 0;
    const int s1Len = profile->s1Len;
    const parasail_matrix_t *matrix = profile->matrix;
    const int32_t segWidth = 2; /* number of values in vector unit */
    const int32_t segLen = (s1Len + segWidth - 1) / segWidth;
    __m128i* const restrict vProfile = (__m128i*)profile->profile64.score;
    __m128i* restrict pvHStore = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHLoad =  parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvE = parasail_memalign___m128i(16, segLen);
    __m128i* restrict pvHMax = parasail_memalign___m128i(16, segLen);
    __m128i vGapO = _mm_set1_epi64x(open);
    __m128i vGapE = _mm_set1_epi64x(gap);
    __m128i vZero = _mm_setzero_si128();
    int64_t score = NEG_INF;
    __m128i vMaxH = vZero;
    __m128i vMaxHUnit = vZero;
    int64_t maxp = INT64_MAX - (int64_t)(matrix->max+1);
    /*int64_t stop = profile->stop == INT32_MAX ?  INT64_MAX : (int64_t)profile->stop;*/
#ifdef PARASAIL_TABLE
    parasail_result_t *result = parasail_result_new_table1(segLen*segWidth, s2Len);
#else
#ifdef PARASAIL_ROWCOL
    parasail_result_t *result = parasail_result_new_rowcol1(segLen*segWidth, s2Len);
    const int32_t offset = (s1Len - 1) % segLen;
    const int32_t position = (segWidth - 1) - (s1Len - 1) / segLen;
#else
#ifdef PARASAIL_TRACE
    __m128i* const restrict pvHT = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvET = parasail_memalign___m128i(16, segLen);
    __m128i* const restrict pvEa = parasail_memalign___m128i(16, segLen);
    parasail_result_t *result = parasail_result_new_trace(segLen*segWidth, s2Len);
    __m128i vTZero = _mm_set1_epi64x(PARASAIL_ZERO);
    __m128i vTIns  = _mm_set1_epi64x(PARASAIL_INS);
    __m128i vTDel  = _mm_set1_epi64x(PARASAIL_DEL);
    __m128i vTDiag = _mm_set1_epi64x(PARASAIL_DIAG);
#else
    parasail_result_t *result = parasail_result_new();
#endif
#endif
#endif

    parasail_memset___m128i(pvHStore, vZero, segLen);
    parasail_memset___m128i(pvE, _mm_set1_epi64x(-open), segLen);

#ifdef PARASAIL_TRACE
    parasail_memset___m128i(pvEa, _mm_set1_epi64x(-open), segLen);
    for (i=0; i<segLen; ++i) {
        arr_store(result->trace_ins_table,
                vTDiag, i, segLen, 0, s2Len);
    }
#endif

    /* outer loop over database sequence */
    for (j=0; j<s2Len; ++j) {
        __m128i vEF_opn;
        __m128i vE;
        __m128i vE_ext;
        __m128i vF;
        __m128i vF_ext;
        __m128i vH;
        __m128i vH_dag;
        const __m128i* vP = NULL;

        /* Initialize F value to 0.  Any errors to vH values will be
         * corrected in the Lazy_F loop. */
        vF = _mm_sub_epi64(vZero,vGapO);

        /* load final segment of pvHStore and shift left by 8 bytes */
        vH = _mm_load_si128(&pvHStore[segLen - 1]);
        vH = _mm_slli_si128(vH, 8);

        /* Correct part of the vProfile */
        vP = vProfile + matrix->mapper[(unsigned char)s2[j]] * segLen;

        if (end_ref == j-2) {
            /* Swap in the max buffer. */
            SWAP3(pvHMax,  pvHLoad,  pvHStore)
        }
        else {
            /* Swap the 2 H buffers. */
            SWAP(pvHLoad,  pvHStore)
        }

        /* inner loop to process the query sequence */
        for (i=0; i<segLen; ++i) {
            vE = _mm_load_si128(pvE + i);

            /* Get max from vH, vE and vF. */
            vH_dag = _mm_add_epi64(vH, _mm_load_si128(vP + i));
            vH_dag = _mm_max_epi64_rpl(vH_dag, vZero);
            vH = _mm_max_epi64_rpl(vH_dag, vE);
            vH = _mm_max_epi64_rpl(vH, vF);
            /* Save vH values. */
            _mm_store_si128(pvHStore + i, vH);

#ifdef PARASAIL_TRACE
            {
                __m128i cond_zero = _mm_cmpeq_epi64(vH, vZero);
                __m128i case1 = _mm_cmpeq_epi64(vH, vH_dag);
                __m128i case2 = _mm_cmpeq_epi64(vH, vF);
                __m128i vT = _mm_blendv_epi8(
                        _mm_blendv_epi8(vTIns, vTDel, case2),
                        _mm_blendv_epi8(vTDiag, vTZero, cond_zero),
                        case1);
                _mm_store_si128(pvHT + i, vT);
                arr_store(result->trace_table, vT, i, segLen, j, s2Len);
            }
#endif
#ifdef PARASAIL_TABLE
            arr_store(result->score_table, vH, i, segLen, j, s2Len);
#endif
            vMaxH = _mm_max_epi64_rpl(vH, vMaxH);
            vEF_opn = _mm_sub_epi64(vH, vGapO);

            /* Update vE value. */
            vE_ext = _mm_sub_epi64(vE, vGapE);
            vE = _mm_max_epi64_rpl(vEF_opn, vE_ext);
            _mm_store_si128(pvE + i, vE);
#ifdef PARASAIL_TRACE
            {
                __m128i vEa = _mm_load_si128(pvEa + i);
                __m128i vEa_ext = _mm_sub_epi64(vEa, vGapE);
                vEa = _mm_max_epi64_rpl(vEF_opn, vEa_ext);
                _mm_store_si128(pvEa + i, vEa);
                if (j+1<s2Len) {
                    __m128i cond = _mm_cmpgt_epi64_rpl(vEF_opn, vEa_ext);
                    __m128i vT = _mm_blendv_epi8(vTIns, vTDiag, cond);
                    _mm_store_si128(pvET + i, vT);
                    arr_store(result->trace_ins_table, vT, i, segLen, j+1, s2Len);
                }
            }
#endif

            /* Update vF value. */
            vF_ext = _mm_sub_epi64(vF, vGapE);
            vF = _mm_max_epi64_rpl(vEF_opn, vF_ext);
#ifdef PARASAIL_TRACE
            {
                __m128i cond = _mm_cmpgt_epi64_rpl(vEF_opn, vF_ext);
                __m128i vT = _mm_blendv_epi8(vTDel, vTDiag, cond);
                if (i+1<segLen) {
                    arr_store(result->trace_del_table, vT, i+1, segLen, j, s2Len);
                }
            }
#endif

            /* Load the next vH. */
            vH = _mm_load_si128(pvHLoad + i);
        }

        /* Lazy_F loop: has been revised to disallow adjecent insertion and
         * then deletion, so don't update E(i, i), learn from SWPS3 */
        for (k=0; k<segWidth; ++k) {
#ifdef PARASAIL_TRACE
            __m128i vFa;
            __m128i vFa_ext;
            __m128i vHp = _mm_load_si128(&pvHLoad[segLen - 1]);
            vHp = _mm_slli_si128(vHp, 8);
            vEF_opn = _mm_slli_si128(vEF_opn, 8);
            vEF_opn = _mm_insert_epi64(vEF_opn, -open, 0);
            vF_ext = _mm_slli_si128(vF_ext, 8);
            vF_ext = _mm_insert_epi64(vF_ext, NEG_INF, 0);
            vFa_ext = vF_ext;
#endif
            vF = _mm_slli_si128(vF, 8);
            vF = _mm_insert_epi64(vF, -open, 0);
#ifdef PARASAIL_TRACE
            vFa = vF;
#endif
            for (i=0; i<segLen; ++i) {
                vH = _mm_load_si128(pvHStore + i);
                vH = _mm_max_epi64_rpl(vH,vF);
                _mm_store_si128(pvHStore + i, vH);
#ifdef PARASAIL_TRACE
                {
                    __m128i vT;
                    __m128i case1;
                    __m128i case2;
                    __m128i cond;
                    vHp = _mm_add_epi64(vHp, _mm_load_si128(vP + i));
                    vHp = _mm_max_epi64_rpl(vHp, vZero);
                    case1 = _mm_cmpeq_epi64(vH, vHp);
                    case2 = _mm_cmpeq_epi64(vH, vF);
                    cond = _mm_andnot_si128(case1,case2);
                    vT = _mm_load_si128(pvHT + i);
                    vT = _mm_blendv_epi8(vT, vTDel, cond);
                    _mm_store_si128(pvHT + i, vT);
                    arr_store(result->trace_table, vT, i, segLen, j, s2Len);
                }
#endif
#ifdef PARASAIL_TABLE
                arr_store(result->score_table, vH, i, segLen, j, s2Len);
#endif
                vMaxH = _mm_max_epi64_rpl(vH, vMaxH);
                /* Update vF value. */
#ifdef PARASAIL_TRACE
                {
                    __m128i cond = _mm_cmpgt_epi64_rpl(vEF_opn, vFa_ext);
                    __m128i vT = _mm_blendv_epi8(vTDel, vTDiag, cond);
                    arr_store(result->trace_del_table, vT, i, segLen, j, s2Len);
                }
#endif
                vEF_opn = _mm_sub_epi64(vH, vGapO);
                vF_ext = _mm_sub_epi64(vF, vGapE);
#ifdef PARASAIL_TRACE
                {
                    __m128i vET = _mm_load_si128(pvET + i);
                    __m128i vEa = _mm_load_si128(pvEa + i);
                    __m128i cond = _mm_cmpgt_epi64_rpl(vEF_opn, vEa);
                    vEa = _mm_max_epi64_rpl(vEa, vEF_opn);
                    _mm_store_si128(pvEa + i, vEa);
                    vET = _mm_blendv_epi8(vET, vTDiag, cond);
                    if (j+1<s2Len) {
                        arr_store(result->trace_ins_table, vET, i, segLen, j+1, s2Len);
                    }
                }
#endif
                if (! _mm_movemask_epi8(
                            _mm_or_si128(
                                _mm_cmpgt_epi64_rpl(vF_ext, vEF_opn),
                                _mm_cmpeq_epi64(vF_ext, vEF_opn))))
                    goto end;
                /*vF = _mm_max_epi64_rpl(vEF_opn, vF_ext);*/
                vF = vF_ext;
#ifdef PARASAIL_TRACE
                vFa_ext = _mm_sub_epi64(vFa, vGapE);
                vFa = _mm_max_epi64_rpl(vEF_opn, vFa_ext);
                vHp = _mm_load_si128(pvHLoad + i);
#endif
            }
        }
end:
        {
        }

#ifdef PARASAIL_ROWCOL
        /* extract last value from the column */
        {
            vH = _mm_load_si128(pvHStore + offset);
            for (k=0; k<position; ++k) {
                vH = _mm_slli_si128(vH, 8);
            }
            result->score_row[j] = (int64_t) _mm_extract_epi64 (vH, 1);
        }
#endif

        {
            __m128i vCompare = _mm_cmpgt_epi64_rpl(vMaxH, vMaxHUnit);
            if (_mm_movemask_epi8(vCompare)) {
                score = _mm_hmax_epi64_rpl(vMaxH);
                /* if score has potential to overflow, abort early */
                if (score > maxp) {
                    result->saturated = 1;
                    break;
                }
                vMaxHUnit = _mm_set1_epi64x(score);
                end_ref = j;
            }
        }

        /*if (score == stop) break;*/
    }

#ifdef PARASAIL_ROWCOL
    for (i=0; i<segLen; ++i) {
        __m128i vH = _mm_load_si128(pvHStore+i);
        arr_store_col(result->score_col, vH, i, segLen);
    }
#endif

    if (score == INT64_MAX) {
        result->saturated = 1;
    }

    if (result->saturated) {
        score = 0;
        end_query = 0;
        end_ref = 0;
    }
    else {
        if (end_ref == j-1) {
            /* end_ref was the last store column */
            SWAP(pvHMax,  pvHStore)
        }
        else if (end_ref == j-2) {
            /* end_ref was the last load column */
            SWAP(pvHMax,  pvHLoad)
        }
        /* Trace the alignment ending position on read. */
        {
            int64_t *t = (int64_t*)pvHMax;
            int32_t column_len = segLen * segWidth;
            end_query = s1Len - 1;
            for (i = 0; i<column_len; ++i, ++t) {
                if (*t == score) {
                    int32_t temp = i / segWidth + i % segWidth * segLen;
                    if (temp < end_query) {
                        end_query = temp;
                    }
                }
            }
        }
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

#ifdef PARASAIL_TRACE
    parasail_free(pvEa);
    parasail_free(pvET);
    parasail_free(pvHT);
#endif
    parasail_free(pvHMax);
    parasail_free(pvE);
    parasail_free(pvHLoad);
    parasail_free(pvHStore);

    return result;
}


