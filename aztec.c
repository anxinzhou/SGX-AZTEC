//
// Created by anxin on 2020/5/9.
//

#include "aztec.h"
#include <stdlib.h>
#include <openssl/evp.h>

G1 h;
G1 ** base_notes_mu;
int ** base_notes_value;
const int y = 1<<26;

void init_aztec_parameters(int max_k_power,int note_ratio){
    if(max_k_power%note_ratio!=0) exit(-2);
    if(1<<max_k_power >= y) exit(-3);
    // init h
    G1_init_field(&h);
    mpz_t rd;
    mpz_init_set_ui(rd,myrand_int());
    G1_mul(&h,&G1Gen,rd);
    mpz_clear(rd);

    // init mu and base note value
    base_notes_mu = (G1**) malloc(sizeof(G1*)*note_ratio);
    base_notes_value = (int **)malloc(sizeof(int*)*note_ratio);
    for(int i=0;i<note_ratio;i++) {
        int base_note_num = max_k_power/note_ratio;

        //BigInteger tmp(y-note_value);
        //tmp = modinv(tmp,BN256::p.getMagnitude());
        //base_note[j] = h * tmp;


        int *notes_value = (int *)malloc(sizeof(int)*(1<<base_note_num));
        notes_value[0] = 0;
        int base_len = 1;
        int k=0;
        while(k<base_note_num) {
            int len = base_len;
            int value = 1<<(i*base_note_num+k);
            for(int j=0;j<len;j++) {
                notes_value[base_len+j] = value + notes_value[j];
            }
            base_len*=2;
            k++;
        }

        G1 *mus = (G1 *)malloc(sizeof(G1)*(1<<base_note_num));
        for(int j=0;j<(1<<base_note_num);j++) {
            mpz_t tmp;
            mpz_init_set_ui(tmp,y-j);
//            tmp = modinv(tmp,BN256::p.getMagnitude());
            //auto t1 = std::chrono::high_resolution_clock::now();
            G1 tmp2;
            G1_init_field(&tmp2);
            G1_mul(&tmp2,&h,tmp);
            //auto t2 = std::chrono::high_resolution_clock::now();
            //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

            //std::cout<<"mul" << duration / 1000.0 / 1000 << "s" << endl;
            mus[j]=tmp2;
            mpz_clear(tmp);
        }
        base_notes_mu[i] = mus;
        base_notes_value[i] = notes_value;
    }
}

Cmt gen_cmt(int k, mpz_t a, G1 *mu) {
    Cmt cmt;
    G1_init_field(&cmt.gamma);
    G1_init_field(&cmt.sigma);

    mpz_t v;
    mpz_init_set_ui(v,k);
    G1_mul(&cmt.gamma,mu,v);

    mpz_t va;
    mpz_init(va);
    mpz_mul(va,v,a);
    G1_mul(&cmt.sigma, mu, va);

    G1 ah;
    G1_init_field(&ah);
    G1_mul(&ah,&h,a);

    G1_add(&cmt.sigma,&cmt.sigma,&ah);

    G1_clear_field(&ah);
    mpz_clears(va,v,NULL);
    return cmt;
}

void cmt_clear_field(Cmt *cmt) {
    G1_clear_field(&cmt->gamma);
    G1_clear_field(&cmt->sigma);
}

Proof gen_proof(Cmt *cmt, int n, int k_public, int *k, mpz_t *a) {
    // check balance
    int balance = 0;
    for(int i=0;i<n;i++) {
        balance += k[i];
    }
    if(balance!=k_public) {
        exit(2);
    }

    // pick ba1, ba2,...,
    // pick bk2,...,bkn
    mpz_t *ba = (mpz_t *)malloc(sizeof(mpz_t)*n);
    mpz_t *bk = (mpz_t *)malloc(sizeof(mpz_t)*n);

    for(int i=0;i<n;i++) {
        mpz_init_set_ui(ba[i],myrand_int());
        if(i==0) {
            mpz_init_set_ui(bk[0],0);
        } else {
            mpz_init_set_ui(bk[i],myrand_int());
            mpz_add(bk[0],bk[0],bk[i]);
        }
    }
    mpz_neg(bk[0],bk[0]);

    // calculate B1 ... Bn
    G1* B = (G1 *)malloc(sizeof(G1)*n);
    G1 tmp;
    G1_init_field(&tmp);
    for(int i=0;i<n;i++) {
        G1_init_field(&B[i]);
        G1_mul(&B[i],&cmt[i].gamma,bk[0]);
        G1_mul(&tmp,&h,ba[i]);
        G1_add(&B[i],&B[i],&tmp);
    }
    G1_clear_field(&tmp);

    // calculate challange
    mpz_t c;
    mpz_init(c);
    gen_challenge(c, cmt, n, B);
    // calculate k1_, ... kn_ and a1_ ... an_
    mpz_t *k_= (mpz_t *)malloc(sizeof(mpz_t)*(n-1));
    mpz_t *a_ = (mpz_t *)malloc(sizeof(mpz_t)*n);
    for(int i=0; i<n; i++) {
        mpz_init_set(a_[i],c);
        mpz_mul(a_[i],a_[i],a[i]);
        mpz_add(a_[i],a_[i],ba[i]);
        if(i==0) continue;
        mpz_init_set(k_[i-1],c);
        mpz_mul_ui(k_[i-1],k_[i-1],k[i]);
        mpz_add(k_[i-1],k_[i-1],bk[i]);
    }

    // free memory
    for(int i=0;i<n;i++) {
        G1_clear_field(&B[i]);
    }
    free(B);
    for(int i=0;i<n;i++) {
        mpz_clears(ba[i],bk[i],NULL);
    }
    free(ba);
    free(bk);

    Proof proof = {.k=k_,.a=a_,.len=n};
    mpz_init_set(proof.c,c);
    mpz_clear(c);
    return proof;
}

void proof_clear_field(Proof *proof) {
    for(int i=0;i<proof->len;i++) {
        if(i==0) mpz_clear(proof->a[i]);
        else
        mpz_clears(proof->a[i],proof->k[i-1],NULL);
    }
    free(proof->a);
    free(proof->k);
}

void sha256(unsigned char *digest, unsigned char *message, size_t message_len) {
    unsigned int SHALEN = 32;
    EVP_MD_CTX *mdctx;
    mdctx = EVP_MD_CTX_create();
    EVP_DigestInit_ex(mdctx, EVP_sha256(), NULL);
    EVP_DigestUpdate(mdctx, message, message_len);
    EVP_DigestFinal_ex(mdctx, digest, &SHALEN);
    EVP_MD_CTX_destroy(mdctx);

//    auto result = ethash::keccak256(message, message_len);
//    digest = result.bytes;
}

void gen_challenge(mpz_t c, Cmt *cmt, int n, G1 *B) {

    int size_of_G1 = 32*2;
    // sizeof(cmt + n + B)
    int message_size = size_of_G1 *2 *n + 32 + size_of_G1 *n;

    unsigned char *message = (unsigned char *)malloc(sizeof(unsigned char)*message_size);
    int start = 0; // record memcpy location

    // copy commitment
    for (int i = 0; i < n; ++i) {
        G1_marshal(message + start, &cmt[i].gamma);
        G1_marshal(message + start + size_of_G1, &cmt[i].sigma);
        start += 2 * size_of_G1;
    }

    // copy n
    // use 32 bytes to be compatible with ethereum uint256
    char padded[32];
    memset(padded,'0',32);
    memcpy(message + start, (unsigned char *) padded, 32 - sizeof(n));
    start += 32 - sizeof(n);
    memcpy(message + start, (unsigned char *) &n, sizeof(n));
    // reverse m to use big end  .... because of little end in linux
    for(int i=0;i<sizeof(n)/2;i++) {
        unsigned char tmp = message[start+i];
        message[start+i] = message[start+sizeof(n)-1-i];
        message[start+sizeof(n)-1-i] = tmp;
    }
    start += sizeof(n);

    //copy B
    for(int i=0; i<n; i++) {
        G1_marshal(message + start, &B[i]);
        start += size_of_G1;
    }


    // hash message
    unsigned char digest[32];
    sha256(digest, message, message_size);

    // convert hash digest to mpz
    char buf[65];
    const char * hex = "0123456789ABCDEF";
    for(int i=0;i<32;i++) {
        uint8_t c = digest[i];
        uint8_t first = c/16;
        uint8_t second = c%16;
        buf[i*2] = hex[(c>>4)&0xf];
        buf[i*2+1] = hex[c&0xf];
    }
    buf[64]='\0';
    mpz_set_str(c,buf,16);
    // clean memory
    free(message);
}