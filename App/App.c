/*
 * Copyright (C) 2011-2020 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

// App.cpp : Define the entry point for the console application.
//
#include <string.h>
#include <assert.h>
//#include <thread>

#include "Enclave_u.h"
#include <sgx_urts.h>
#include "sgx_tseal.h"
#include <stdio.h>

//#include <omp.h>

#define ENCLAVE_NAME "enclave.signed.so"

// Global data
sgx_enclave_id_t global_eid = 0;


typedef struct _sgx_errlist_t {
	sgx_status_t err;
	const char *msg;
	const char *sug; /* Suggestion */
} sgx_errlist_t;


/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] =
		{ { SGX_ERROR_UNEXPECTED, "Unexpected error occurred.",
		NULL }, { SGX_ERROR_INVALID_PARAMETER, "Invalid parameter.",
		NULL }, { SGX_ERROR_OUT_OF_MEMORY, "Out of memory.",
		NULL }, { SGX_ERROR_ENCLAVE_LOST, "Power transition occurred.",
				"Please refer to the sample \"PowerTransition\" for details." },
				{ SGX_ERROR_INVALID_ENCLAVE, "Invalid enclave image.",
				NULL }, { SGX_ERROR_INVALID_ENCLAVE_ID,
						"Invalid enclave identification.",
						NULL }, { SGX_ERROR_INVALID_SIGNATURE,
						"Invalid enclave signature.",
						NULL }, { SGX_ERROR_OUT_OF_EPC, "Out of EPC memory.",
				NULL },
				{ SGX_ERROR_NO_DEVICE, "Invalid SGX device.",
						"Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards." },
				{ SGX_ERROR_MEMORY_MAP_CONFLICT, "Memory map conflicted.",
				NULL }, { SGX_ERROR_INVALID_METADATA,
						"Invalid enclave metadata.",
						NULL }, { SGX_ERROR_DEVICE_BUSY, "SGX device was busy.",
				NULL }, { SGX_ERROR_INVALID_VERSION,
						"Enclave version was invalid.",
						NULL }, { SGX_ERROR_INVALID_ATTRIBUTE,
						"Enclave was not authorized.",
						NULL }, { SGX_ERROR_ENCLAVE_FILE_ACCESS,
						"Can't open enclave file.",
						NULL }, };

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret) {
	size_t idx = 0;
	size_t ttl = sizeof sgx_errlist / sizeof sgx_errlist[0];

	for (idx = 0; idx < ttl; idx++) {
		if (ret == sgx_errlist[idx].err) {
			if (NULL != sgx_errlist[idx].sug)
				printf("Info: %s\n", sgx_errlist[idx].sug);
			printf("Error: %s\n", sgx_errlist[idx].msg);
			break;
		}
	}

	if (idx == ttl)
		printf(
				"Error code is 0x%X. Please refer to the \"Intel SGX SDK Developer Reference\" for more details.\n",
				ret);
}

/* OCall functions */
void ocall_print_string(const char *str) {
	/* Proxy/Bridge will check the length and null-terminate
	 * the input string to prevent buffer overflow.
	 */
	printf("%s", str);
}



void ocall_get_time(unsigned long int* t){
	*t= clock();
}

int main(int argc, char *argv[]) {
	//omp_set_dynamic(0);     // Explicitly disable dynamic teams
	//omp_set_num_threads(8); // Use 4 threads for all consecutive
	(void) argc, (void) argv;


	sgx_status_t ret = sgx_create_enclave(ENCLAVE_NAME, 1,
			NULL, NULL, &global_eid,
			NULL);
	if (ret != SGX_SUCCESS) {
		print_error_message(ret);
		return -1;
	}
/*
	// test k power impact
	int k_sample = 9;
	int user_num = 4000;
	int note_ratio= 2;
	int counts_of_k[9] = {16,18,20,22,24,26,28,30,32};


	printf("%d\n",sizeof(unsigned long int));
	for(int i=0;i<k_sample;i++) {

		double ob_time,nonob_time;
		obnote_generation_benchmark(global_eid,counts_of_k[i],user_num,note_ratio,&ob_time);
		nobnote_generation_benchmark(global_eid,counts_of_k[i],user_num,&nonob_time);
		printf("ratio:%d kpower:%d user:%d ob time %f\n",note_ratio,counts_of_k[i],user_num,ob_time);
		printf("ratio:%d kpower:%d user:%d non ob time %f\n",note_ratio,counts_of_k[i],user_num,nonob_time);
	}
*/

	// test note ratio impact
int note_sample = 11;
	int max_k_power = 32;
	int user_num = 4000;
	int counts_of_ratio[11] = {2,3,4,5,6,7,8,10,12,14,16};
	//int counts_of_ratio[11] = {6,7,8,10,12,14,16};

	printf("%d\n",sizeof(unsigned long int));
	for(int i=0;i<note_sample;i++) {

		double ob_time,nonob_time;
		obnote_generation_benchmark(global_eid,max_k_power,user_num,counts_of_ratio[i],&ob_time);
		nobnote_generation_benchmark(global_eid,max_k_power,user_num,&nonob_time);
		printf("ratio:%d kpower:%d user:%d ob time %f\n",counts_of_ratio[i],max_k_power,user_num,ob_time);
		printf("ratio:%d kpower:%d user:%d non ob time %f\n",counts_of_ratio[i],max_k_power,user_num,nonob_time);
	}



	// test user num
	/*int user_sample = 10;
	int counts_of_user[user_sample];
	for (int i=0;i<user_sample;i++) {
		counts_of_user[i] = (i+1)*1000;
	}

	printf("%d\n",sizeof(unsigned long int));
	for(int i=0;i<user_sample;i++) {
		int note_ratio = 2;
		int max_k_power = 28;
		//int user_num = 1000;
		double ob_time,nonob_time;
		obnote_generation_benchmark(global_eid,max_k_power,counts_of_user[i],note_ratio,&ob_time);
		nobnote_generation_benchmark(global_eid,max_k_power,counts_of_user[i],&nonob_time);
		printf("ratio:%d kpower:%d user:%d ob time %f\n",note_ratio,max_k_power,counts_of_user[i],ob_time);
		printf("ratio:%d kpower:%d user:%d non ob time %f\n",note_ratio,max_k_power,counts_of_user[i],nonob_time);
	}*/

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

	//test(global_eid);
	sgx_destroy_enclave(global_eid);

	//getchar();
	return 0;
}

