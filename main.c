//
// Created by anxin on 2020/5/6.
//

#include "bn256.h"
#include <stdio.h>
#include <time.h>
#include "math.h"
#include "util.h"
#include "aztec.h"
#include <stdlib.h>



// max_k is 10, 12
void obnote_generation_benchmark(int max_k_power,int user_num,int note_ratio) {
    if(max_k_power%note_ratio!=0) {
        printf("wrong note_ratio\n");
        exit(-4);
    }

    int *user_original_reward = (int *)malloc(sizeof(int)*user_num);
    int k_public = 0;
    for(int i=0; i<user_num; i++) {
        user_original_reward[i] = myrand_int()%(1<<max_k_power);
        k_public += user_original_reward[i];
    }
    mpz_t *user_a = (mpz_t *)malloc(sizeof(mpz_t)*user_num*note_ratio);
    for(int i=0;i<user_num*note_ratio;i++) {
        mpz_init_set_ui(user_a[i],myrand_int());
    }

    // benchmark start from here
    clock_t t1 = clock();
    //critical part, get mu  for each user
    G1 *user_mus = (G1 *) malloc(sizeof(G1)*user_num*note_ratio);
    int *user_values = (int *)malloc(sizeof(int)*user_num*note_ratio);
    for(int i=0;i<user_num;i++) {
        int reward = user_original_reward[i];
        unsigned int binary_ind = reward;

//        // calculate binary
//        for(int j=max_k_power-1;j>=0;j--) {
//            uint8_t flag = reward >= (1<<j);
//            binary_ind |= (1<<j) * flag;
//            reward -= (1<<j)*flag;
//        }
        // calcualte the sum in each division
        int * division_value = (int*) malloc(sizeof(int)*note_ratio);
        for(int j=0;j<note_ratio;j++) {
            int base_num = max_k_power / note_ratio;
            int total = 0;
            for(int k=0;k< base_num;k++) {
                int pos = j*base_num + k;
                unsigned int value = 1<<pos;
                total+= value * ( (binary_ind >> (pos))&1);
            }
            division_value[j] = total;
        }
        // linear scan each division
        for(int j=0; j<note_ratio; j++) {
            int note_value = division_value[j];
            int base_num = max_k_power / note_ratio;
            G1 mu;
//            G1 test_mu;
            G1_init_field(&mu);
//            G1_init_field(&test_mu);
//            bool find = false; //test
            for(int k=0;k<(1<<base_num);++k) {
                G1 * base_mu = &base_notes_mu[j][k];
                int base_value = base_notes_value[j][k];
                oblivious_assign_G1(base_value==note_value,&mu,base_mu,&mu);
//                if(base_value==note_value) { //test
//                    find = true;
//                    G1_set_g1(&test_mu,&base_mu);
//                }
            }
//            if(!find) throw("not found");
//            if(mpz_cmp(mu.x,test_mu.x)!=0) {
//                cout<<"wrong mu"<<endl;
//                return;
//            }
            //if(!(mu+(-test_mu)).isInfinity()) throw("oblivious assingment is not right");
            user_mus[i*note_ratio+j] = mu;
            user_values[i*note_ratio+j]= note_value;
        }
        free(division_value);
    }
    // prepare proof
    int note_num = user_num * note_ratio;
    Cmt *user_cmt = (Cmt *)malloc(sizeof(Cmt)*note_num);
    for(int i=0; i<note_num; i++) {
        Cmt cmt = gen_cmt(user_values[i],user_a[i],&user_mus[i]);
        user_cmt[i] = cmt;
    }
    Proof proof = gen_proof(user_cmt,note_num,k_public,user_values,user_a);
    clock_t t2 = clock();
    double duration = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("proof generation time %f \n",duration);

    // clear memory
    proof_clear_field(&proof);
    for(int i=0; i <note_num;i++) {
        cmt_clear_field(&user_cmt[i]);
    }
    free(user_cmt);
    for(int i=0;i<note_num;i++) {
        mpz_clear(user_a[i]);
    }
    free(user_a);
    for(int i=0; i< note_num;i++) {
        G1_clear_field(&user_mus[i]);
    }
    free(user_values);
    free(user_mus);
    free(user_original_reward);
}

void nobnote_generation_benchmark(int max_k_power,int user_num) {


    int k_public = 0;
    mpz_t *user_a = (mpz_t *)malloc(sizeof(mpz_t)*user_num);
    for(int i=0;i<user_num;i++) {
        mpz_init_set_ui(user_a[i],myrand_int());
    }

    // benchmark start from here
    clock_t t1 = clock();
    //critical part, get mu  for each user
    G1 *user_mus = (G1 *) malloc(sizeof(G1)*user_num);
    int *user_values = (int *)malloc(sizeof(int)*user_num);
    for(int i=0;i<user_num;i++) {
        user_values[i] = myrand_int()%(1<<max_k_power);
        k_public += user_values[i];

        mpz_t tmp;
        mpz_init_set_ui(tmp,y-user_values[i]);
        G1 mu;
        G1_init_field(&mu);
        G1_mul(&mu,&h,tmp);

        user_mus[i] = mu;
        mpz_clear(tmp);
    }
    // prepare proof
    int note_num = user_num;
    Cmt *user_cmt = (Cmt *)malloc(sizeof(Cmt)*note_num);
    for(int i=0; i<note_num; i++) {
        Cmt cmt = gen_cmt(user_values[i],user_a[i],&user_mus[i]);
        user_cmt[i] = cmt;
    }
    Proof proof = gen_proof(user_cmt,note_num,k_public,user_values,user_a);
    clock_t t2 = clock();
    double duration = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("proof generation time %f \n",duration);

    // clear memory
    proof_clear_field(&proof);
    for(int i=0; i <note_num;i++) {
        cmt_clear_field(&user_cmt[i]);
    }
    free(user_cmt);
    for(int i=0;i<note_num;i++) {
        mpz_clear(user_a[i]);
    }
    free(user_a);
    for(int i=0; i< note_num;i++) {
        G1_clear_field(&user_mus[i]);
    }
    free(user_values);
    free(user_mus);
}

void g1_test() {
    init_curve_parameter();
    G1 a;
    G1_init_field(&a);
    mpz_t rd;
    mpz_init_set_ui(rd,myrand_int());
    G1_set_g1(&a,&G1Gen);
    G1 b;
    G1_init_field(&b);
    G1_set_g1(&b,&a);
    mpz_t tmp;
    mpz_init_set_ui(tmp,4);

    G1_mul(&a,&a,tmp);
//    G1_add(&a,&a,&G1Gen);
//    G1_add(&a,&a,&G1Gen);
//    G1_add(&a,&a,&G1Gen);

    G1_add(&b,&b,&b);
    G1_add(&b,&b,&b);
    G1_negate(&b);


    G1_add(&a,&a,&b);
    printf("%d\n",G1_is_infinity(&a));
    //obnote_generation_benchmark(8,1,2);
}

void g2_test() {
    init_curve_parameter();
    G2 a;
    G2_init_field(&a);
    mpz_t rd;
    mpz_init_set_ui(rd,myrand_int());
    G2_set_g2(&a,&G2Gen);
    G2 b;
    G2_init_field(&b);
    G2_set_g2(&b,&a);
    mpz_t tmp;
    mpz_init_set_ui(tmp,8);

//    G2_mul(&a,&a,tmp);
    G2_add(&a,&a,&G2Gen);
    G2_add(&a,&a,&G2Gen);
    G2_add(&a,&a,&G2Gen);
//
    G2_add(&b,&b,&b);
    G2_add(&b,&b,&b);
    G2_negate(&b);


    G2_add(&a,&a,&b);
    printf("%d\n",G2_is_infinity(&a));
    //obnote_generation_benchmark(8,1,2);
}

int main() {
    int note_ratio = 2;
    int max_k_power = 20;
    int user_num = 100;
    init_curve_parameter();
    init_aztec_parameters(max_k_power,note_ratio);
    obnote_generation_benchmark(max_k_power,user_num,note_ratio);
    nobnote_generation_benchmark(max_k_power,user_num);
//    mpz_t t;
//    mpz_init_set_str(t,"ffff",16);
//    mpz_t t2;
//    mpz_init_set_str(t2,"-000ffff",16);
//    char st[64];
//    char st2[64];
//    mpz_get_str(st,16,t);
//    mpz_get_str(st2,16,t2);
//    cout<<st<<endl;
//    cout<<st2<<endl;
}
