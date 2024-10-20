#include <stdio.h>

#include "common.h"


void test(int a, float* floates) {
    printf("Hello : %d\n", a);
}

__device__ const int CAL_SAMPLES = 426;
__global__ void fit(float *X, float *y, float *XTest) {

    __shared__ float sharedY[CAL_SAMPLES];
    __shared__ float sharedWeight;
    __shared__ float diff[CAL_SAMPLES];


    __shared__ int epochs;
    __shared__ int currentEpoch;
    __shared__ float learningRate;
    __shared__ float bias;

    int tid = threadIdx.x * blockDim.y + threadIdx.y;

    if(tid < CAL_SAMPLES) {
        sharedY[tid] = y[tid];
    }

    if(threadIdx.x == 0 && threadIdx.y == 0) {
        learningRate = 0.001;
        bias = 1.0;
        epochs = 1000;
        sharedWeight = 0.1;
        currentEpoch = 0;
    }

    __syncthreads();

    while(currentEpoch < epochs) {

        if(tid < CAL_SAMPLES) {
            float prediction = 1.0 / (1.0 + exp(-(X[tid] * sharedWeight + bias)));

            diff[tid] = (prediction - sharedY[tid]);

        }

        __syncthreads();

        if(threadIdx.x == 0 && threadIdx.y == 0) {
            float sum = 0.0;
            float diffSum = 0.0;
            for(int i = 0; i < CAL_SAMPLES; i++) {
                sum += diff[i] * X[i];
                diffSum += diff[i];
            }

            sharedWeight -= learningRate * ((1.0 / (float) CAL_SAMPLES) * sum);
            //printf("BRUH : %f\n", learningRate * (1.0 / (float) CAL_SAMPLES * sum));
            bias -= learningRate * ((1.0 / (float) CAL_SAMPLES) * diffSum);
            currentEpoch++;
        }

        __syncthreads();
    }


    if(tid < 114) {
        float prediction = 1.0 / (1.0 + exp(-(XTest[tid] * sharedWeight + bias)));

        printf("%d[%f] = %f\n", tid, XTest[tid], prediction);

    }

    if(tid == 0) {
        printf("Weight : %f - Bias : %f\n", sharedWeight, bias);
    }

}

__global__ void fitCalibrationModel(float *X, int *y, int* mask, float *weights, float *globalBias) {

    __shared__ int sharedY[SAMPLES];
    __shared__ float diff[SAMPLES];
    __shared__ float sharedWeight;
    __shared__ float bias;
    __shared__ float learningRate;
    __shared__ int currentEpoch;
    __shared__ int epochs;

    int tid = threadIdx.x * blockDim.y + threadIdx.y;
    printf("%d\n", tid);
    if(tid < SAMPLES) {
        int index = tid;
        sharedY[index] = y[mask[index]] == blockIdx.x;
        index += (blockDim.x * blockDim.y);

        sharedY[index] = y[mask[index]] == blockIdx.x;
        index += (blockDim.x * blockDim.y);

        if(index < SAMPLES ) {
            sharedY[index] = y[mask[index]] == blockIdx.x;
        }
    }

    if(threadIdx.x == 0 && threadIdx.y == 0) {
        epochs = 1000;
        currentEpoch = 0;
        learningRate = 0.01;
        bias = globalBias[blockIdx.x];
        sharedWeight = weights[blockIdx.x];
    }

    __syncthreads();

    while(currentEpoch < epochs) {
        tid = threadIdx.x * blockDim.y + threadIdx.y;
        for(int k = 0; k < 3; k++) {
            if (tid < SAMPLES) {
                float prediction = 1.0 / (1.0 + exp(-(X[MAX_CLASS * tid + blockIdx.x] * sharedWeight + bias)));
                diff[tid] = (prediction - (float) sharedY[tid]);
            }

            tid += blockDim.x * blockDim.y;
        }
        __syncthreads();


        if(threadIdx.x == 0 && threadIdx.y == 0) {
            float sum = 0.0;
            float diffSum = 0.0;
            for(int i = 0; i < SAMPLES; i++) {
                sum += diff[i] * X[MAX_CLASS * i + blockIdx.x];
                diffSum += diff[i];
            }

            sharedWeight -= learningRate * ((1.0 / (float) SAMPLES) * sum);
            bias -= learningRate * ((1.0 / (float) SAMPLES) * diffSum);
            currentEpoch++;
        }

        __syncthreads();
    }

    if(tid == 0) {
        weights[blockIdx.x] = sharedWeight;
        globalBias[blockIdx.x] = bias;

        printf("Final weight : %f ----- Final bias : %f\n", sharedWeight, bias);
    }

}


void cudaFit(float *X, float *y, float *test, int sampleSize, int testSize, int featureSize) {


    float* d_x;
    float* d_x_test;
    float* d_y;

    cudaMalloc(&d_x, sizeof(float) * sampleSize * featureSize);
    cudaMemcpy(d_x, X, sizeof(float) * sampleSize * featureSize, cudaMemcpyHostToDevice);

    cudaMalloc(&d_y, sizeof(float) * sampleSize);
    cudaMemcpy(d_y, y, sizeof(float) * sampleSize, cudaMemcpyHostToDevice);

    cudaMalloc(&d_x_test, sizeof(float) * testSize * featureSize);
    cudaMemcpy(d_x_test, test, sizeof(float) * testSize * featureSize, cudaMemcpyHostToDevice);

    dim3 gridSize(1,1);
    dim3 blockSize(sampleSize,1);

    fit<<<gridSize, blockSize>>>(d_x, d_y, d_x_test);
    cudaDeviceSynchronize();
}

