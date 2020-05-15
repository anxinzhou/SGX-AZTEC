//
// Created by anxin on 2020/5/9.
//

#ifndef PROOF2_AZTEC_H
#define PROOF2_AZTEC_H

#include "bn256.h"
#include "util.h"

extern G1 h;
extern G1 ** base_notes_mu;
extern int ** base_notes_value;
extern unsigned long int y;

void ob_init_aztec_parameters(int max_k_power,int note_ratio);
void non_ob_init_aztec_parameters();

void non_ob_clear_aztec_parameters();

void ob_clear_aztec_parameters(int max_k_power, int note_ratio);

struct Commitment {
    G1 gamma;
    G1 sigma;
};

typedef struct Commitment Cmt;

Cmt gen_cmt(int k, mpz_t a, G1 *mu);
void cmt_clear_field(Cmt *cmt);

struct Proof {
    mpz_t *a;
    mpz_t *k;
    mpz_t c;
    int len;
};

typedef struct Proof Proof;

Proof gen_proof(Cmt *cmt, int n, int k_public, int *k, mpz_t *a);

void proof_clear_field(Proof *proof);

void gen_challenge(mpz_t c, Cmt *cmt, int n, G1 *B);
#endif //PROOF2_AZTEC_H
