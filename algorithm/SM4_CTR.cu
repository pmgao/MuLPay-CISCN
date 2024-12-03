#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <chrono>

// Device constant S-box for non-linear transformation in SM4 encryption on GPU
__device__ __constant__ uint8_t d_Sbox[256] = {
	0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05,
	0x2b, 0x67, 0x9a, 0x76, 0x2a, 0xbe, 0x04, 0xc3, 0xaa, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
	0x9c, 0x42, 0x50, 0xf4, 0x91, 0xef, 0x98, 0x7a, 0x33, 0x54, 0x0b, 0x43, 0xed, 0xcf, 0xac, 0x62,
	0xe4, 0xb3, 0x1c, 0xa9, 0xc9, 0x08, 0xe8, 0x95, 0x80, 0xdf, 0x94, 0xfa, 0x75, 0x8f, 0x3f, 0xa6,
	0x47, 0x07, 0xa7, 0xfc, 0xf3, 0x73, 0x17, 0xba, 0x83, 0x59, 0x3c, 0x19, 0xe6, 0x85, 0x4f, 0xa8,
	0x68, 0x6b, 0x81, 0xb2, 0x71, 0x64, 0xda, 0x8b, 0xf8, 0xeb, 0x0f, 0x4b, 0x70, 0x56, 0x9d, 0x35,
	0x1e, 0x24, 0x0e, 0x5e, 0x63, 0x58, 0xd1, 0xa2, 0x25, 0x22, 0x7c, 0x3b, 0x01, 0x21, 0x78, 0x87,
	0xd4, 0x00, 0x46, 0x57, 0x9f, 0xd3, 0x27, 0x52, 0x4c, 0x36, 0x02, 0xe7, 0xa0, 0xc4, 0xc8, 0x9e,
	0xea, 0xbf, 0x8a, 0xd2, 0x40, 0xc7, 0x38, 0xb5, 0xa3, 0xf7, 0xf2, 0xce, 0xf9, 0x61, 0x15, 0xa1,
	0xe0, 0xae, 0x5d, 0xa4, 0x9b, 0x34, 0x1a, 0x55, 0xad, 0x93, 0x32, 0x30, 0xf5, 0x8c, 0xb1, 0xe3,
	0x1d, 0xf6, 0xe2, 0x2e, 0x82, 0x66, 0xca, 0x60, 0xc0, 0x29, 0x23, 0xab, 0x0d, 0x53, 0x4e, 0x6f,
	0xd5, 0xdb, 0x37, 0x45, 0xde, 0xfd, 0x8e, 0x2f, 0x03, 0xff, 0x6a, 0x72, 0x6d, 0x6c, 0x5b, 0x51,
	0x8d, 0x1b, 0xaf, 0x92, 0xbb, 0xdd, 0xbc, 0x7f, 0x11, 0xd9, 0x5c, 0x41, 0x1f, 0x10, 0x5a, 0xd8,
	0x0a, 0xc1, 0x31, 0x88, 0xa5, 0xcd, 0x7b, 0xbd, 0x2d, 0x74, 0xd0, 0x12, 0xb8, 0xe5, 0xb4, 0xb0,
	0x89, 0x69, 0x97, 0x4a, 0x0c, 0x96, 0x77, 0x7e, 0x65, 0xb9, 0xf1, 0x09, 0xc5, 0x6e, 0xc6, 0x84,
	0x18, 0xf0, 0x7d, 0xec, 0x3a, 0xdc, 0x4d, 0x20, 0x79, 0xee, 0x5f, 0x3e, 0xd7, 0xcb, 0x39, 0x48
};

// Host-side S-box used for key expansion
const uint8_t h_Sbox[256] = {
	0xd6, 0x90, 0xe9, 0xfe, 0xcc, 0xe1, 0x3d, 0xb7, 0x16, 0xb6, 0x14, 0xc2, 0x28, 0xfb, 0x2c, 0x05,
	0x2b, 0x67, 0x9a, 0x76, 0x2a, 0xbe, 0x04, 0xc3, 0xaa, 0x44, 0x13, 0x26, 0x49, 0x86, 0x06, 0x99,
	0x9c, 0x42, 0x50, 0xf4, 0x91, 0xef, 0x98, 0x7a, 0x33, 0x54, 0x0b, 0x43, 0xed, 0xcf, 0xac, 0x62,
	0xe4, 0xb3, 0x1c, 0xa9, 0xc9, 0x08, 0xe8, 0x95, 0x80, 0xdf, 0x94, 0xfa, 0x75, 0x8f, 0x3f, 0xa6,
	0x47, 0x07, 0xa7, 0xfc, 0xf3, 0x73, 0x17, 0xba, 0x83, 0x59, 0x3c, 0x19, 0xe6, 0x85, 0x4f, 0xa8,
	0x68, 0x6b, 0x81, 0xb2, 0x71, 0x64, 0xda, 0x8b, 0xf8, 0xeb, 0x0f, 0x4b, 0x70, 0x56, 0x9d, 0x35,
	0x1e, 0x24, 0x0e, 0x5e, 0x63, 0x58, 0xd1, 0xa2, 0x25, 0x22, 0x7c, 0x3b, 0x01, 0x21, 0x78, 0x87,
	0xd4, 0x00, 0x46, 0x57, 0x9f, 0xd3, 0x27, 0x52, 0x4c, 0x36, 0x02, 0xe7, 0xa0, 0xc4, 0xc8, 0x9e,
	0xea, 0xbf, 0x8a, 0xd2, 0x40, 0xc7, 0x38, 0xb5, 0xa3, 0xf7, 0xf2, 0xce, 0xf9, 0x61, 0x15, 0xa1,
	0xe0, 0xae, 0x5d, 0xa4, 0x9b, 0x34, 0x1a, 0x55, 0xad, 0x93, 0x32, 0x30, 0xf5, 0x8c, 0xb1, 0xe3,
	0x1d, 0xf6, 0xe2, 0x2e, 0x82, 0x66, 0xca, 0x60, 0xc0, 0x29, 0x23, 0xab, 0x0d, 0x53, 0x4e, 0x6f,
	0xd5, 0xdb, 0x37, 0x45, 0xde, 0xfd, 0x8e, 0x2f, 0x03, 0xff, 0x6a, 0x72, 0x6d, 0x6c, 0x5b, 0x51,
	0x8d, 0x1b, 0xaf, 0x92, 0xbb, 0xdd, 0xbc, 0x7f, 0x11, 0xd9, 0x5c, 0x41, 0x1f, 0x10, 0x5a, 0xd8,
	0x0a, 0xc1, 0x31, 0x88, 0xa5, 0xcd, 0x7b, 0xbd, 0x2d, 0x74, 0xd0, 0x12, 0xb8, 0xe5, 0xb4, 0xb0,
	0x89, 0x69, 0x97, 0x4a, 0x0c, 0x96, 0x77, 0x7e, 0x65, 0xb9, 0xf1, 0x09, 0xc5, 0x6e, 0xc6, 0x84,
	0x18, 0xf0, 0x7d, 0xec, 0x3a, 0xdc, 0x4d, 0x20, 0x79, 0xee, 0x5f, 0x3e, 0xd7, 0xcb, 0x39, 0x48
};

// System parameters (CK) for SM4 algorithm
const uint32_t ck[32] = {
	0x00070e15, 0x1c232a31, 0x383f464d, 0x545b6269,
	0x70777e85, 0x8c939aa1, 0xa8afb6bd, 0xc4cbd2d9,
	0xe0e7eef5, 0xfc030a11, 0x181f262d, 0x343b4249,
	0x50575e65, 0x6c737a81, 0x888f969d, 0xa4abb2b9,
	0xc0c7ced5, 0xdce3eaf1, 0xf8ff060d, 0x141b2229,
	0x30373e45, 0x4c535a61, 0x686f767d, 0x848b9299,
	0xa0a7aeb5, 0xbcc3cad1, 0xd8dfe6ed, 0xf4fb0209,
	0x10171e25, 0x2c333a41, 0x484f565d, 0x646b7279 };

// Basic bitwise operation macros
#define Rotl(_x, _y) (((_x) << (_y)) | ((_x) >> (32 - (_y))))  
#define L1(_B) ((_B) ^ Rotl(_B, 2) ^ Rotl(_B, 10) ^ Rotl(_B, 18) ^ Rotl(_B, 24))  
#define L2(_B) ((_B) ^ Rotl(_B, 13) ^ Rotl(_B, 23))

// Host-side S-box transformation macro
#define h_S_box(_A) (h_Sbox[(_A) >> 24 & 0xFF] << 24 ^ \
					h_Sbox[(_A) >> 16 & 0xFF] << 16 ^ \
					h_Sbox[(_A) >>  8 & 0xFF] <<  8 ^ \
					h_Sbox[(_A) & 0xFF])  
#define d_S_box(_A) (d_Sbox[(_A) >> 24 & 0xFF] << 24 ^ \
					d_Sbox[(_A) >> 16 & 0xFF] << 16 ^ \
					d_Sbox[(_A) >>  8 & 0xFF] <<  8 ^ \
					d_Sbox[(_A) & 0xFF])  

// Core encryption function running on GPU device
// Processes a single 16-byte block using SM4 algorithm in CTR mode
__device__ __inline__ void SM4Crypt(const uint8_t* input, uint8_t* output, uint32_t* rk, const uint8_t* plaintext) {
    uint32_t mid, x0, x1, x2, x3;

	// Convert input bytes to 32-bit words with big-endian ordering
	x0 = (input[3] << 24) | (input[2] << 16) | (input[1] << 8) | (input[0] << 0); 
	x1 = (input[7] << 24) | (input[6] << 16) | (input[5] << 8) | (input[4] << 0); 
	x2 = (input[11] << 24) | (input[10] << 16) | (input[9] << 8) | (input[8] << 0); 
	x3 = (input[15] << 24) | (input[14] << 16) | (input[13] << 8) | (input[12] << 0); 

	// 32 rounds of SM4 encryption
    #pragma unroll
    for (uint32_t r = 0; r < 32; r += 4) {
        mid = x1 ^ x2 ^ x3 ^ rk[r + 0];
        mid = d_S_box(mid);
        x0 ^= L1(mid);
        mid = x2 ^ x3 ^ x0 ^ rk[r + 1];
        mid = d_S_box(mid);
        x1 ^= L1(mid);
        mid = x3 ^ x0 ^ x1 ^ rk[r + 2];
        mid = d_S_box(mid);
        x2 ^= L1(mid);
        mid = x0 ^ x1 ^ x2 ^ rk[r + 3];
        mid = d_S_box(mid);
        x3 ^= L1(mid);
    }

	// Convert 32-bit words back to bytes with reverse ordering for output
    output[0] = (x3 >> 24) & 0xFF;
    output[1] = (x3 >> 16) & 0xFF;
    output[2] = (x3 >> 8) & 0xFF;
    output[3] = (x3 >> 0) & 0xFF;

    output[4] = (x2 >> 24) & 0xFF;
    output[5] = (x2 >> 16) & 0xFF;
    output[6] = (x2 >> 8) & 0xFF;
    output[7] = (x2 >> 0) & 0xFF;

    output[8] = (x1 >> 24) & 0xFF;
    output[9] = (x1 >> 16) & 0xFF;
    output[10] = (x1 >> 8) & 0xFF;
    output[11] = (x1 >> 0) & 0xFF;

    output[12] = (x0 >> 24) & 0xFF;
    output[13] = (x0 >> 16) & 0xFF;
    output[14] = (x0 >> 8) & 0xFF;
    output[15] = (x0 >> 0) & 0xFF;

	// XOR with plaintext for CTR mode operation
    #pragma unroll
    for (int i = 0; i < 16; i++) {
        output[i] ^= plaintext[i];
    }
}

// CUDA kernel for SM4 encryption in CTR mode
// Processes multiple blocks in parallel
__global__ void SM4_encrypt_cuda_ctr(uint8_t* p, uint8_t* c, uint32_t* rk, int num_cipher_blocks, uint8_t* counter) {
	int i = (blockIdx.x * blockDim.x) + threadIdx.x;
	if (i < num_cipher_blocks) {
		uint8_t ca[16];

		// Copy counter value and increment for this block
		#pragma unroll
		for (int w = 1; w < 16; w++) {
			ca[w] = counter[w];
		}
		ca[15] += i; 
		
		SM4Crypt(ca, &(c[i * 16]), rk, &(p[i * 16]));
	}
}

// Key expansion function for SM4
// Generates round keys from the initial key
void SM4KeyExt(uint8_t* Key, uint32_t* rk) {
	uint32_t r, mid, x0, x1, x2, x3;

	// Convert key bytes to 32-bit words
	x0 = (Key[3] << 24) | (Key[2] << 16) | (Key[1] << 8) | (Key[0] << 0); 
	x1 = (Key[7] << 24) | (Key[6] << 16) | (Key[5] << 8) | (Key[4] << 0); 
	x2 = (Key[11] << 24) | (Key[10] << 16) | (Key[9] << 8) | (Key[8] << 0); 
	x3 = (Key[15] << 24) | (Key[14] << 16) | (Key[13] << 8) | (Key[12] << 0); 

	// XOR with system parameters
	x0 ^= 0xa3b1bac6;
	x1 ^= 0x56aa3350;
	x2 ^= 0x677d9197;
	x3 ^= 0xb27022dc;

	// Generate 32 round keys
	#pragma unroll
	for (r = 0; r < 32; r += 4) {
		mid = x1 ^ x2 ^ x3 ^ ck[r];
		mid = h_S_box(mid);
		rk[r] = x0 ^= L2(mid);
		mid = x2 ^ x3 ^ x0 ^ ck[r + 1];
		mid = h_S_box(mid);
		rk[r + 1] = x1 ^= L2(mid);
		mid = x3 ^ x0 ^ x1 ^ ck[r + 2];
		mid = h_S_box(mid);
		rk[r + 2] = x2 ^= L2(mid);
		mid = x0 ^ x1 ^ x2 ^ ck[r + 3];
		mid = h_S_box(mid);
		rk[r + 3] = x3 ^= L2(mid);
	}
}

// Main CUDA wrapper function for SM4 encryption
// Handles device memory allocation, data transfer, and kernel launch
cudaError_t SM4_cuda_ctr(char* text, char* cipher, int size, uint32_t* key, uint8_t* Counter) {
	cudaError_t cudaStatus;
	uint8_t* p, * c, * counter;
	uint32_t* rk;
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, 0);

	int blocks, threads;
	int n = size >> 4;

	// Select GPU device
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cuda SetDevice failed\n");
		goto Error;
	}

	// Allocate device memory for counter, plaintext, ciphertext, and round keys
	cudaStatus = cudaMalloc((void**)&counter, 16 * sizeof(uint8_t));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&p, size * sizeof(uint8_t));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&c, size * sizeof(uint8_t));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed\n");
		goto Error;
	}
	cudaStatus = cudaMalloc((void**)&rk, 32 * sizeof(uint32_t));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed\n");
		goto Error;
	}

	// Copy data from host to device
	cudaStatus = cudaMemcpy(counter, Counter, 16 * sizeof(uint8_t), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed");
		goto Error;
	}

	cudaStatus = cudaMemcpy(p, text, size * sizeof(uint8_t), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed");
		goto Error;
	}

	cudaStatus = cudaMemcpy(rk, key, 32 * sizeof(uint32_t), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed");
		goto Error;
	}

	// Calculate grid and block dimensions for kernel launch
	if (n > 1) {
		threads = (n < prop.maxThreadsPerBlock) ? n : prop.maxThreadsPerBlock;
		blocks = n / threads;
		if (n % threads != 0) {
			blocks += 1;
		}
	}
	else {
		threads = 1; blocks = 1;
	}

	// launch kernel
	SM4_encrypt_cuda_ctr <<< blocks, threads >>> (p, c, rk, n, counter);

	// Check for kernel launch errors
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Kernel launch failed:%s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// Wait for kernel completion
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching kernel\n", cudaStatus);
		goto Error;
	}

	// Copy results back to host
	cudaStatus = cudaMemcpy(cipher, c, size * sizeof(uint8_t), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed\n");
		goto Error;
	}


Error:
// Free device memory
	cudaFree(c);
	cudaFree(p);
	cudaFree(rk);
	return cudaStatus;
}

int main() {
	// Initialize counter for CTR mode
	uint8_t Counter[16] = { 0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10 };

	// Read encryption key from file
	FILE* file_;
	char* filename;
	unsigned char key[16];
	filename = (char*)"key.txt";
	file_ = fopen(filename, "r");
	if (file_ == NULL) {
		printf("open file failed!\n");
		exit(1);
	}
	for (size_t i = 0; i < 16; i++) {
		char temp[3];
		temp[0] = fgetc(file_);
		temp[1] = fgetc(file_);
		temp[3] = 0;
		key[i] = std::stoi(temp, nullptr, 16);
	}
	fclose(file_);

	// Generate round keys
	uint32_t rk[32];
	SM4KeyExt(key, rk);

	// Read input file
	std::ifstream inFILE("message.txt", std::ios::binary);
	if (!inFILE) {
		printf("Failed to open the file\n");
		exit(1);
	}
	int L = inFILE.tellg();
	inFILE.seekg(0, std::ios::end);
	int M = inFILE.tellg();
	int size = M-L;
	inFILE.seekg(0, inFILE.beg);

	// Add padding if necessary
	int add;
	add = (16 - size % 16) % 16;
	char* p = (char*)calloc(size + add, sizeof(char));
	memset(p, 0, size + add);
	
	// Read input data
	inFILE.read(p, size);
	inFILE.close();
	size = size + add;
	
	// launch GPU
	char* m = (char*)calloc(size, sizeof(char));
	SM4_cuda_ctr(p, m, size, rk, Counter);

	// Write output data
	std::ofstream FILE("message.txt", std::ios::binary);
	if (!FILE) {
		printf("Failed to open the file\n");
		exit(1);
	}
	FILE.write(m, size);
	FILE.close();

	// free host memory
	free(p);
	free(m);
	return 0;
}