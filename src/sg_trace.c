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

#include "parasail.h"
#include "parasail/memory.h"

#define NEG_INF_32 (INT32_MIN/2)
#define MAX(a,b) ((a)>(b)?(a):(b))

#define ENAME parasail_sg_trace

parasail_result_t* ENAME(
        const char * const restrict _s1, const int s1Len,
        const char * const restrict _s2, const int s2Len,
        const int open, const int gap, const parasail_matrix_t *matrix)
{
    parasail_result_t *result = parasail_result_new_trace(s1Len, s2Len);
    int * const restrict s1 = parasail_memalign_int(16, s1Len);
    int * const restrict s2 = parasail_memalign_int(16, s2Len);
    int * const restrict tbl_pr = parasail_memalign_int(16, s2Len+1);
    int * const restrict del_pr = parasail_memalign_int(16, s2Len+1);
    int i = 0;
    int j = 0;
    int score = NEG_INF_32;
    int end_query = s1Len;
    int end_ref = s2Len;

    for (i=0; i<s1Len; ++i) {
        s1[i] = matrix->mapper[(unsigned char)_s1[i]];
    }
    for (j=0; j<s2Len; ++j) {
        s2[j] = matrix->mapper[(unsigned char)_s2[j]];
    }

    /* upper left corner */
    tbl_pr[0] = 0;
    del_pr[0] = NEG_INF_32;
    
    /* first row */
    for (j=1; j<=s2Len; ++j) {
        tbl_pr[j] = 0;
        del_pr[j] = NEG_INF_32;
    }

    /* iter over first sequence */
    for (i=1; i<s1Len; ++i) {
        const int * const restrict matrow = &matrix->matrix[matrix->size*s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        for (j=1; j<=s2Len; ++j) {
            int del_tbl;
            int del_del;
            int ins_tbl;
            int ins_ins;
            int tbl_tbl;
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_tbl = Nscore - open;
            del_del = del_pr[j] - gap;
            ins_tbl = Wscore - open;
            ins_ins = ins_cr    - gap;
            tbl_tbl = NWscore + matrow[s2[j-1]];
            del_pr[j] = MAX(del_tbl, del_del);
            ins_cr    = MAX(ins_tbl, ins_ins);
            Wscore = MAX(tbl_tbl, del_pr[j]);
            Wscore = MAX(Wscore, ins_cr);
            tbl_pr[j] = Wscore;
            result->trace_del_table[(i-1)*s2Len + (j-1)] = 
                (del_tbl > del_del) ? PARASAIL_DIAG
                                    : PARASAIL_DEL;
            result->trace_ins_table[(i-1)*s2Len + (j-1)] = 
                (ins_tbl > ins_ins) ? PARASAIL_DIAG
                                    : PARASAIL_INS;
            result->trace_table[(i-1)*s2Len + (j-1)] = 
                (Wscore == tbl_tbl) ? PARASAIL_DIAG
                    : (Wscore == del_pr[j]) ? PARASAIL_DEL
                    : PARASAIL_INS;
        }
        if (Wscore > score) {
            score = Wscore;
            end_query = i-1;
            end_ref = s2Len-1;
        }
    }
    {
        /* i == s1Len */
        const int * const restrict matrow = &matrix->matrix[matrix->size*s1[i-1]];
        /* init first column */
        int Nscore = tbl_pr[0];
        int Wscore = 0;
        int ins_cr = NEG_INF_32;
        tbl_pr[0] = Wscore;
        for (j=1; j<=s2Len; ++j) {
            int del_tbl;
            int del_del;
            int ins_tbl;
            int ins_ins;
            int tbl_tbl;
            int NWscore = Nscore;
            Nscore = tbl_pr[j];
            del_tbl = Nscore - open;
            del_del = del_pr[j] - gap;
            ins_tbl = Wscore - open;
            ins_ins = ins_cr    - gap;
            tbl_tbl = NWscore + matrow[s2[j-1]];
            del_pr[j] = MAX(del_tbl, del_del);
            ins_cr    = MAX(ins_tbl, ins_ins);
            Wscore = MAX(tbl_tbl, del_pr[j]);
            Wscore = MAX(Wscore, ins_cr);
            tbl_pr[j] = Wscore;
            if (Wscore > score) {
                score = Wscore;
                end_query = s1Len-1;
                end_ref = j-1;
            }
            else if (Wscore == score && j-1 < end_ref) {
                end_query = s1Len-1;
                end_ref = j-1;
            }
            result->trace_del_table[(i-1)*s2Len + (j-1)] = 
                (del_tbl > del_del) ? PARASAIL_DIAG
                                    : PARASAIL_DEL;
            result->trace_ins_table[(i-1)*s2Len + (j-1)] = 
                (ins_tbl > ins_ins) ? PARASAIL_DIAG
                                    : PARASAIL_INS;
            result->trace_table[(i-1)*s2Len + (j-1)] = 
                (Wscore == tbl_tbl) ? PARASAIL_DIAG
                    : (Wscore == del_pr[j]) ? PARASAIL_DEL
                    : PARASAIL_INS;
        }
    }

    result->score = score;
    result->end_query = end_query;
    result->end_ref = end_ref;

    parasail_free(del_pr);
    parasail_free(tbl_pr);
    parasail_free(s2);
    parasail_free(s1);

    return result;
}

