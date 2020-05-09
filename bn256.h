//
// Created by anxin on 2020/5/5.
//

#ifndef PROOF_BN256_H
#define PROOF_BN256_H

#include <gmp.h>
#include "util.h"
#include <string.h>
#include <stdint.h>

struct G1;
struct G2;
struct gfP2;
typedef struct G1 G1;
typedef struct G2 G2;
typedef struct gfP2 gfP2;
extern mpz_t q;
extern G1 G1Gen;
extern G2 G2Gen;
extern const int y;

void init_curve_parameter();

struct G1 {
    mpz_t x;
    mpz_t y;
    mpz_t z;
    mpz_t t;
};

struct gfP2 {
    mpz_t x;
    mpz_t y;
};

struct G2 {
    gfP2 x;
    gfP2 y;
    gfP2 z;
    gfP2 t;
};


//void G1_init(G1 **p);
//void G1_clear(G1 **p);
void G1_init_field(G1 *p);

void G1_clear_field(G1 *p);

void G1_add(G1 *c, G1 *a, G1 *b);

void G1_mul(G1 *c, G1 *a, mpz_t b);

int G1_is_infinity(G1 *p);

void G1_set_infinity(G1 *p);

void G1_double(G1 *p);

void G1_negate(G1 *p);

void G1_set_g1(G1 *p, G1 *a);

void G1_make_affine(G1 *p);

void G1_marshal(unsigned char *buf, G1 *p);

//void G2_init(G2 **p);
//void G2_clear(G2 **p);
void G2_init_field(G2 *p);

void G2_clear_field(G2 *p);

void G2_set_g2(G2 *p, G2 *a);

void G2_add(G2 *c, G2 *a, G2 *b);

void G2_mul(G2 *c, G2 *a, mpz_t b);

int G2_is_infinity(G2 *p);

void G2_set_infinity(G2 *p);

void G2_double(G2 *p);

void G2_negate(G2 *p);

//void gfP2_init(gfP2 ** p);
//void gfP2_init_field(gfP2 **p);
void gfP2_clear(gfP2 *p);

void gfP2_clear_field(gfP2 *p);

void gfP2_add(gfP2 *c, gfP2 *a, gfP2 *b);

void gfP2_sub(gfP2 *c, gfP2 *a, gfP2 *b);

void gfP2_mul(gfP2 *c, gfP2 *a, gfP2 *b);

int gfP2_is_zero(gfP2 *p);

void gfP2_set_zero(gfP2 *p);

int gfP2_is_one(gfP2 *p);

void gfP2_set_one(gfP2 *p);

void gfP2_invert(gfP2 *p);

void gfP2_square(gfP2 *p);

void gfP2_set_gfP2(gfP2 *p, gfP2 *a);


void oblivious_assign_str(uint8_t flag, char *res, char *t_val, char *f_val, int len);

void oblivious_assign_mpz(uint8_t flag, mpz_t res, mpz_t t_val, mpz_t f_val);

void oblivious_assign_G1(uint8_t flag, G1 *res, G1 *t_val, G1 *f_val);

#endif //PROOF_BN256_H
