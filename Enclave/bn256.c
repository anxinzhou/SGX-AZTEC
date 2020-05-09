//
// Created by anxin on 2020/5/5.
//

#include "bn256.h"
#include <stdlib.h>

mpz_t q;
G1 G1Gen;
G2 G2Gen;

void init_curve_parameter() {
    // init g1
    mpz_init_set_ui(G1Gen.x, 1);
    mpz_init_set_si(G1Gen.y, -2);
    mpz_init_set_ui(G1Gen.z, 1);
    mpz_init_set_ui(G1Gen.t, 1);

    // init g2gen
    mpz_init_set_str(G2Gen.x.x, "21167961636542580255011770066570541300993051739349375019639421053990175267184", 10);
    mpz_init_set_str(G2Gen.x.y, "64746500191241794695844075326670126197795977525365406531717464316923369116492", 10);

    mpz_init_set_str(G2Gen.y.x, "20666913350058776956210519119118544732556678129809273996262322366050359951122", 10);
    mpz_init_set_str(G2Gen.y.y, "17778617556404439934652658462602675281523610326338642107814333856843981424549", 10);

    mpz_init_set_str(G2Gen.z.x, "0", 10);
    mpz_init_set_str(G2Gen.z.y, "1", 10);

    mpz_init_set_str(G2Gen.t.x, "0", 10);
    mpz_init_set_str(G2Gen.t.y, "1", 10);


    mpz_init_set_str(q, "65000549695646603732796438742359905742825358107623003571877145026864184071783", 10);
}

void clear_curve_parameter() {
	G1_clear_field(&G1Gen);
	G2_clear_field(&G2Gen);
	mpz_clear(q);
}

// init field of p
void G1_init(G1 *p) {
    p = (G1 *) malloc(sizeof(G1));
    mpz_inits(p->x, p->y, p->z, p->t, NULL);
}

void G1_init_field(G1 *p) {
    mpz_inits(p->x, p->y, p->z, p->t, NULL);
}

void G1_clear(G1 *p) {
    mpz_clears(p->x, p->y, p->z, p->t, NULL);
    free(p);
}

void G1_clear_field(G1 *p) {
    mpz_clears(p->x, p->y, p->z, p->t, NULL);
}

// set p to a;
void G1_set_g1(G1 *p, G1 *a) {
    if (p == a) return;
    mpz_set(p->x, a->x);
    mpz_set(p->y, a->y);
    mpz_set(p->z, a->z);
    mpz_set(p->t, a->t);
}


int G1_is_infinity(G1 *p) {
    return mpz_sgn(p->z) == 0;
}

void G1_set_infinity(G1 *p) {
    mpz_set_ui(p->x, 0);
    mpz_set_ui(p->y, 0);
    mpz_set_ui(p->z, 0);
    mpz_set_ui(p->t, 0);
}

void G1_add(G1 *c, G1 *a, G1 *b) {
    if (G1_is_infinity(a)) {
        G1_set_g1(c, b);
        return;
    }
    if (G1_is_infinity(b)) {
        G1_set_g1(c, a);
        return;
    }

    mpz_t z1z1, z2z2, u1, u2, t, s1, s2, h, i, j, r, v, t4, t6;
    mpz_inits(z1z1, z2z2, u1, u2, t, s1, s2, h, i, j, r, v, t4, t6, NULL);

    mpz_mul(z1z1, a->z, a->z);
    mpz_mod(z1z1, z1z1, q);
    mpz_mul(z2z2, b->z, b->z);
    mpz_mod(z2z2, z2z2, q);
    mpz_mul(u1, a->x, z2z2);
    mpz_mod(u1, u1, q);
    mpz_mul(u2, b->x, z1z1);
    mpz_mod(u2, u2, q);

    mpz_mul(t, b->z, z2z2);
    mpz_mod(t, t, q);
    mpz_mul(s1, a->y, t);
    mpz_mod(s1, s1, q);

    mpz_mul(t, a->z, z1z1);
    mpz_mod(t, t, q);
    mpz_mul(s2, b->y, t);
    mpz_mod(s2, s2, q);

    mpz_sub(h, u2, u1);
    int xEqual = mpz_sgn(h) == 0;
    mpz_add(t, h, h);
    mpz_mul(i, t, t);
    mpz_mod(i, i, q);
    mpz_mul(j, h, i);
    mpz_mod(j, j, q);
    mpz_sub(t, s2, s1);

    int yEqual = mpz_sgn(t) == 0;
    if (xEqual && yEqual) {
        G1_set_g1(c, a);
        G1_double(c);
        return;
    }

    mpz_add(r, t, t);
    mpz_mul(v, u1, i);
    mpz_mod(v, v, q);

    mpz_mul(t4, r, r);
    mpz_mod(t4, t4, q);
    mpz_add(t, v, v);
    mpz_sub(t6, t4, j);
    mpz_sub(c->x, t6, t);

    mpz_sub(t, v, c->x);
    mpz_mul(t4, s1, j);
    mpz_mod(t4, t4, q);
    mpz_add(t6, t4, t4);
    mpz_mul(t4, r, t);
    mpz_mod(t4, t4, q);
    mpz_sub(c->y, t4, t6);

    mpz_add(t, a->z, b->z);
    mpz_mul(t4, t, t);
    mpz_mod(t4, t4, q);
    mpz_sub(t, t4, z1z1);
    mpz_sub(t4, t, z2z2);
    mpz_mul(c->z, t4, h);
    mpz_mod(c->z, c->z, q);

    mpz_clears(z1z1, z2z2, u1, u2, t, s1, s2, h, i, j, r, v, t4, t6, NULL);
}

void G1_mul(G1 *c, G1 *a, mpz_t b) {

    G1 sum, t;
    G1_init_field(&t);
    G1_init_field(&sum);

    G1_set_infinity(&sum);
    int len = mpz_sizeinbase(b, 2);
    for (int i = len; i >= 0; i--) {
        G1_set_g1(&t, &sum);
        G1_double(&t);
        if (mpz_tstbit(b, i) != 0) {
            G1_add(&sum, &t, a);
        } else {
            G1_set_g1(&sum, &t);
        }
    }
    G1_set_g1(c, &sum);
    G1_clear_field(&t);
    G1_clear_field(&sum);
}

void G1_double(G1 *p) {
    mpz_t A, B, C, t, t2, d, e, f;
    mpz_inits(A, B, C, t, t2, d, e, f, NULL);

    mpz_mul(A, p->x, p->x);
    mpz_mod(A, A, q);
    mpz_mul(B, p->y, p->y);
    mpz_mod(B, B, q);
    mpz_mul(C, B, B);
    mpz_mod(C, C, q);

    mpz_add(t, p->x, B);
    mpz_mul(t2, t, t);
    mpz_mod(t2, t2, q);
    mpz_sub(t, t2, A);
    mpz_sub(t2, t, C);
    mpz_add(d, t2, t2);
    mpz_add(t, A, A);
    mpz_add(e, t, A);
    mpz_mul(f, e, e);
    mpz_mod(f, f, q);

    mpz_add(t, d, d);
    mpz_sub(p->x, f, t);

    mpz_add(t, C, C);
    mpz_add(t2, t, t);
    mpz_add(t, t2, t2);


    // keep old y
    mpz_set(A, p->y);

    mpz_sub(p->y, d, p->x);
    mpz_mul(t2, e, p->y);
    mpz_mod(t2, t2, q);
    mpz_sub(p->y, t2, t);

    // add old y here
    mpz_mul(t, A, p->z);
    //mpz_mul(t,p->y,p->z);
    mpz_mod(t, t, q);
    mpz_add(p->z, t, t);
    mpz_clears(A, B, C, t, t2, d, e, f, NULL);
}

void G1_negate(G1 *p) {
    mpz_neg(p->y, p->y);
    mpz_set_ui(p->t, 0);
}

void G1_make_affine(G1 *p) {
    if (mpz_sizeinbase(p->z, 2) == 1 && mpz_tstbit(p->z, 0) == 1) {
        return;
    }

    if (G1_is_infinity(p)) {
        mpz_set_ui(p->x, 0);
        mpz_set_ui(p->y, 1);
        mpz_set_ui(p->z, 0);
        mpz_set_ui(p->t, 0);
        return;
    }

    mpz_t zInv, t, zInv2;
    mpz_inits(zInv, t, zInv2, NULL);
    mpz_invert(zInv, p->z, q);
    mpz_mul(t, p->y, zInv);
    mpz_mod(t, t, q);
    mpz_mul(zInv2, zInv, zInv);
    mpz_mod(zInv2, zInv2, q);
    mpz_mul(p->y, t, zInv2);
    mpz_mod(p->y, p->y, q);
    mpz_mul(t, p->x, zInv2);
    mpz_mod(t, t, q);
    mpz_set(p->x, t);
    mpz_set_ui(p->z, 1);
    mpz_set_ui(p->t, 1);
    mpz_clears(zInv, t, zInv2, NULL);
}



uint64_t xtou64(const char *str, int len)
{
    uint64_t res = 0;
    char c;
    int i=0;
    while (i<len) {
        c = *str++;
        char v = (c & 0xF) + (c >> 6) | ((c >> 3) & 0x8);
        res = (res << 4) | (uint64_t) v;
        i+=1;
    }

    return res;
}


void G1_marshal(unsigned char *buf, G1 *p) {
    //assume buf has at least size of 65

    G1 pcopy;
    G1_init_field(&pcopy);
    G1_set_g1(&pcopy, p);

    G1_make_affine(&pcopy);
    if (mpz_sgn(pcopy.x) < 0) {
        mpz_mod(pcopy.x, pcopy.x, q);
    }
    if (mpz_sgn(pcopy.y) < 0) {
        mpz_mod(pcopy.y, pcopy.y, q);
    }

    // padding with prefix 0
    char tmp_buf[65];
    mpz_get_str(tmp_buf, 16, pcopy.x);
    int len_x = strlen(tmp_buf);
    if (len_x < 64) {
        memmove(tmp_buf + 64 - len_x, tmp_buf, len_x);
        memset(tmp_buf, '0', 64 - len_x);
    }
    tmp_buf[64] = '\0';

    int step = sizeof(uint64_t) * 2;
    for (int i = 0; i < 64 / step; i += 1) {
        //clear x;
        uint64_t x = xtou64(tmp_buf+i*step,step);
        *(uint64_t*)(buf+i*step/2) = x;
    }

    mpz_get_str(tmp_buf, 16, pcopy.y);
    int len_y = strlen(tmp_buf);
    if (len_y < 64) {
        memmove(tmp_buf + 64 - len_y, tmp_buf, len_x);
        memset(tmp_buf, '0', 64 - len_y);
    }
    tmp_buf[64] = '\0';
    int shift = 32;
    for (int i = 0; i < 64 / step; i += 1) {
        //clear x;
        uint64_t x = xtou64(tmp_buf+i*step,step);
        *(uint64_t*)(buf+shift+i*step/2) = x;
    }

    G1_clear_field(&pcopy);
}

void gfP2_init(gfP2 *p) {
    p = (gfP2 *) malloc(sizeof(gfP2));
    mpz_inits(p->x, p->y, NULL);
}

void gfP2_init_field(gfP2 *p) {
    mpz_inits(p->x, p->y, NULL);
}

void gfP2_add(gfP2 *c, gfP2 *a, gfP2 *b) {
    mpz_add(c->x, a->x, b->x);
    mpz_add(c->y, a->y, b->y);
}

void gfP2_sub(gfP2 *c, gfP2 *a, gfP2 *b) {
    mpz_sub(c->x, a->x, b->x);
    mpz_sub(c->y, a->y, b->y);
}

void gfP2_mul(gfP2 *c, gfP2 *a, gfP2 *b) {
    mpz_t tx, t, ty;
    mpz_inits(tx, t, ty, NULL);
    mpz_mul(tx, a->x, b->y);
    mpz_mul(t, b->x, a->y);
    mpz_add(tx, tx, t);
    mpz_mod(tx, tx, q);

    mpz_mul(ty, a->y, b->y);
    mpz_mul(t, a->x, b->x);
    mpz_sub(ty, ty, t);
    mpz_mod(c->y, ty, q);
    mpz_set(c->x, tx);
}

int gfP2_is_zero(gfP2 *p) {
    return mpz_sgn(p->x) == 0 && mpz_sgn(p->y) == 0;
}

void gfP2_set_zero(gfP2 *p) {
    mpz_set_ui(p->x, 0);
    mpz_set_ui(p->y, 0);
}

int gfP2_is_one(gfP2 *p) {
    if (mpz_sgn(p->x) == 0) return 0;
    return mpz_sizeinbase(p->y, 2) == 1 && mpz_tstbit(p->y, 1) == 1;
}

void gfP2_set_one(gfP2 *p) {
    mpz_set_ui(p->x, 0);
    mpz_set_ui(p->y, 1);
}

void gfP2_invert(gfP2 *p) {
    mpz_t t, t2, inv;
    mpz_inits(t, t2, inv, NULL);

    mpz_mul(t, p->y, p->y);
    mpz_mul(t2, p->x, p->x);
    mpz_add(t, t, t2);
    mpz_invert(inv, t, q);

    mpz_neg(p->x, p->x);
    mpz_mul(p->x, p->x, inv);
    mpz_mod(p->x, p->x, q);

    mpz_mul(p->y, p->y, inv);
    mpz_mod(p->y, p->y, q);

    mpz_clears(t, t2, inv, NULL);
}

void gfP2_square(gfP2 *p) {
    mpz_t t1, t2, ty;
    mpz_inits(t1, t2, ty, NULL);
    mpz_sub(t1, p->y, p->x);
    mpz_add(t2, p->x, p->y);
    mpz_mul(ty, t1, t2);
    mpz_mod(ty, ty, q);

    mpz_mul(t1, p->x, p->y);
    mpz_add(t1, t1, t1);

    mpz_mod(p->x, t1, q);
    mpz_set(p->y, ty);
    mpz_clears(t1, t2, ty, NULL);
}

void gfP2_set_gfP2(gfP2 *p, gfP2 *a) {
    if (p == a) return;
    mpz_set(p->x, a->x);
    mpz_set(p->y, a->y);
}

void gfP2_clear(gfP2 *p) {
    mpz_clears(p->x, p->y, NULL);
    free(p);
}

void gfP2_clear_field(gfP2 *p) {
    mpz_clears(p->x, p->y, NULL);
}

void G2_init(G2 *p) {
    p = (G2 *) malloc(sizeof(G2));
    mpz_inits(p->x.x, p->x.y,
              p->y.x, p->y.y,
              p->z.x, p->z.y,
              p->t.x, p->t.y, NULL);
}

void G2_init_field(G2 *p) {
    mpz_inits(p->x.x, p->x.y,
              p->y.x, p->y.y,
              p->z.x, p->z.y,
              p->t.x, p->t.y, NULL);
}

void G2_clear_field(G2 *p) {
    mpz_clears(p->x.x, p->x.y,
               p->y.x, p->y.y,
               p->z.x, p->z.y,
               p->t.x, p->t.y, NULL);
}

int G2_is_infinity(G2 *p) {
    return gfP2_is_zero(&p->z);
}

void G2_set_infinity(G2 *p) {
    gfP2_set_zero(&p->x);
    gfP2_set_zero(&p->y);
    gfP2_set_zero(&p->z);
    gfP2_set_zero(&p->t);
}

void G2_double(G2 *p) {
    gfP2 A, B, C, t, t2, d, e, f;
    gfP2_init_field(&A);
    gfP2_init_field(&B);
    gfP2_init_field(&C);
    gfP2_init_field(&t);
    gfP2_init_field(&t2);
    gfP2_init_field(&d);
    gfP2_init_field(&e);
    gfP2_init_field(&f);

    gfP2_set_gfP2(&A, &p->x);
    gfP2_set_gfP2(&B, &p->y);
    gfP2_square(&A);
    gfP2_square(&B);
    gfP2_set_gfP2(&C, &B);
    gfP2_square(&C);

    gfP2_add(&t, &p->x, &B);
    gfP2_set_gfP2(&t2, &t);
    gfP2_square(&t2);
    gfP2_sub(&t, &t2, &A);
    gfP2_sub(&t2, &t, &C);
    gfP2_add(&d, &t2, &t2);
    gfP2_add(&t, &A, &A);
    gfP2_add(&e, &t, &A);
    gfP2_set_gfP2(&f, &e);
    gfP2_square(&f);

    gfP2_add(&t, &d, &d);
    gfP2_sub(&p->x, &f, &t);

    gfP2_add(&t, &C, &C);
    gfP2_add(&t2, &t, &t);
    gfP2_add(&t, &t2, &t2);

    // keep old y
    gfP2_set_gfP2(&A, &p->y);

    gfP2_sub(&p->y, &d, &p->x);
    gfP2_mul(&t2, &e, &p->y);
    gfP2_sub(&p->y, &t2, &t);

    // use old y
    gfP2_mul(&t, &A, &p->z);
    //gfP2_mul(&t,&p->y,&p->z);
    gfP2_add(&p->z, &t, &t);

    gfP2_clear_field(&A);
    gfP2_clear_field(&B);
    gfP2_clear_field(&C);
    gfP2_clear_field(&t);
    gfP2_clear_field(&t2);
    gfP2_clear_field(&d);
    gfP2_clear_field(&e);
    gfP2_clear_field(&f);
}

void G2_negate(G2 *p) {
    gfP2 tmp;
    gfP2_init_field(&tmp);

    gfP2_set_zero(&tmp);
    gfP2_sub(&p->y, &tmp, &p->y);

    gfP2_set_zero(&p->t);
    gfP2_clear_field(&tmp);
}

void G2_add(G2 *c, G2 *a, G2 *b) {
    if (G2_is_infinity(a)) {
        G2_set_g2(c, b);
        return;
    }
    if (G2_is_infinity(b)) {
        G2_set_g2(c, a);
        return;
    }

    gfP2 z1z1, z2z2, u1, u2, t, s1, s2, h, i, j, r, v, t4, t6;
    gfP2_init_field(&z1z1);
    gfP2_init_field(&z2z2);
    gfP2_init_field(&u1);
    gfP2_init_field(&u2);
    gfP2_init_field(&t);
    gfP2_init_field(&s1);
    gfP2_init_field(&s2);
    gfP2_init_field(&h);
    gfP2_init_field(&i);
    gfP2_init_field(&j);
    gfP2_init_field(&r);
    gfP2_init_field(&v);
    gfP2_init_field(&t4);
    gfP2_init_field(&t6);

    gfP2_set_gfP2(&z1z1, &a->z);
    gfP2_set_gfP2(&z2z2, &b->z);
    gfP2_square(&z1z1);
    gfP2_square(&z2z2);
    gfP2_mul(&u1, &a->x, &z2z2);
    gfP2_mul(&u2, &b->x, &z1z1);

    gfP2_mul(&t, &b->z, &z2z2);
    gfP2_mul(&s1, &a->y, &t);

    gfP2_mul(&t, &a->z, &z1z1);
    gfP2_mul(&s2, &b->y, &t);

    gfP2_sub(&h, &u2, &u1);
    int xEqual = gfP2_is_zero(&h);

    gfP2_add(&t, &h, &h);
    gfP2_set_gfP2(&i, &t);
    gfP2_square(&i);
    gfP2_mul(&j, &h, &i);
    gfP2_sub(&t, &s2, &s1);
    int yEqual = gfP2_is_zero(&t);
    if (xEqual && yEqual) {
        G2_set_g2(c, a);
        G2_double(c);
        return;
    }

    gfP2_add(&r, &t, &t);
    gfP2_mul(&v, &u1, &i);
    gfP2_set_gfP2(&t4, &r);
    gfP2_square(&t4);
    gfP2_add(&t, &v, &v);
    gfP2_sub(&t6, &t4, &j);
    gfP2_sub(&c->x, &t6, &t);

    gfP2_sub(&t, &v, &c->x);
    gfP2_mul(&t4, &s1, &j);
    gfP2_add(&t6, &t4, &t4);
    gfP2_mul(&t4, &r, &t);
    gfP2_sub(&c->y, &t4, &t6);

    gfP2_add(&t, &a->z, &b->z);
    gfP2_set_gfP2(&t4, &t);
    gfP2_square(&t4);
    gfP2_sub(&t, &t4, &z1z1);
    gfP2_sub(&t4, &t, &z2z2);
    gfP2_mul(&c->z, &t4, &h);

    gfP2_clear_field(&z1z1);
    gfP2_clear_field(&z2z2);
    gfP2_clear_field(&u1);
    gfP2_clear_field(&u2);
    gfP2_clear_field(&t);
    gfP2_clear_field(&s1);
    gfP2_clear_field(&s2);
    gfP2_clear_field(&i);
    gfP2_clear_field(&j);
    gfP2_clear_field(&h);
    gfP2_clear_field(&r);
    gfP2_clear_field(&v);
    gfP2_clear_field(&t4);
    gfP2_clear_field(&t6);
}

void G2_mul(G2 *c, G2 *a, mpz_t b) {
    G2 sum, t;
    G2_init_field(&t);
    G2_init_field(&sum);

    G2_set_infinity(&sum);
    int len = mpz_sizeinbase(b, 2);
    for (int i = len; i >= 0; i--) {
        G2_set_g2(&t, &sum);
        G2_double(&t);
        if (mpz_tstbit(b, i) != 0) {
            G2_add(&sum, &t, a);
        } else {
            G2_set_g2(&sum, &t);
        }
    }
    G2_set_g2(c, &sum);
    G2_clear_field(&t);
    G2_clear_field(&sum);
}

void G2_set_g2(G2 *p, G2 *a) {
    if (p == a) return;
    gfP2_set_gfP2(&p->x, &a->x);
    gfP2_set_gfP2(&p->y, &a->y);
    gfP2_set_gfP2(&p->z, &a->z);
    gfP2_set_gfP2(&p->t, &a->t);
}

void G2_Clear(G2 *p) {
    mpz_clears(p->x.x, p->x.y,
               p->y.x, p->y.y,
               p->z.x, p->z.y,
               p->t.x, p->t.y, NULL);
    free(p);
}


void oblivious_assign_str(uint8_t flag, char *res, char *t_val, char *f_val, int len) {
    for (int i = 0; i < len; i += 4) {
        uint32_t v = oblivious_assign_CMOV(flag, *(uint32_t *) (t_val + i), *(uint32_t *) (f_val + i));
        *(uint32_t *) (res + i) = v;
    }
    res[len] = '\0';
}

void oblivious_assign_G1(uint8_t flag, G1 *res, G1 *t_val, G1 *f_val) {
    oblivious_assign_mpz(flag, res->x, t_val->x, f_val->x);
    oblivious_assign_mpz(flag, res->y, t_val->y, f_val->y);
    oblivious_assign_mpz(flag, res->z, t_val->z, f_val->z);
    oblivious_assign_mpz(flag, res->t, t_val->t, f_val->t);
}

void oblivious_assign_mpz(uint8_t flag, mpz_t res, mpz_t t_val, mpz_t f_val) {
    char ts[69]; //1bit sign +  64bits data + 1bit terminal +  extra 3bits for padding
    char fs[69];
    char rs[69]; // extra 3bits for padding
    mpz_get_str(ts, 16, t_val);
    mpz_get_str(fs, 16, f_val);
    // pad to 64bit
    int len_t = strlen(ts);
    int len_f = strlen(fs);
    int max_l = len_t > len_f ? len_t : len_f;
    if (max_l % 4 != 0) max_l = (max_l / 4 + 1) * 4;

    if (len_t < max_l) {
        int shift = 0;
        if (ts[0] == '-') shift = 1;
        memmove(shift + ts + max_l - len_t, shift + ts, len_t);
        memset(ts + shift, '0', max_l - len_t);
    }
    if (len_f < max_l) {
        int shift = 0;
        if (fs[0] == '-') shift = 1;
        memmove(shift + fs + max_l - len_f, shift + fs, len_f);
        memset(shift + fs, '0', max_l - len_f);
    }
    oblivious_assign_str(flag, rs, ts, fs, max_l);
    mpz_set_str(res, rs, 16);
}
