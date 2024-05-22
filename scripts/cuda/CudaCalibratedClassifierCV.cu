#include "../common/common.h"
#include <stdio.h>

#include "CudaCalibratedClassifierCV.h"


__device__ const int MAX_CLASS = 92;
__device__ const int TREE_NUM = 2000;
__device__ const int MAX_DEPTH = 8;

__device__ const int SAMPLES = 2801;
__device__ int FEATURE_SIZE = 1000;
__device__ const int TRAINING_FOLD_SIZE = 2240;
__device__ const int VALIDATION_FOLD_SIZE = 560;
__device__ const int LEAF_NODES_PER_TREE = 256;
__device__ const int FEATURE_SQRT = 96;
__device__ const int MIN_SAMPLES_SPLIT = 2;
__device__ const int MIN_SAMPLES_PER_LEAF = 1;
__device__ const float IG_THRESHOlD = 0.1;


dim3 gridSize(TREE_NUM, 1);
dim3 blockSize(10, FEATURE_SQRT);

dim3 predictionGrid(1, 1);
dim3 predictionBlock(1024, 1);

dim3 calibrationGrid(MAX_CLASS, 1);
dim3 calibrationBlock(VALIDATION_FOLD_SIZE, 1);

cudaError_t cudaStatus;

void checkCudaError(cudaError_t cudaError, std::string functionName) {
    cudaError = cudaGetLastError();

    if (cudaError != cudaSuccess) {
        std::cout << "(" << functionName << ") " << "CUDA error : " << cudaGetErrorString(cudaError) << std::endl;
        exit(1);
    }

}

void printPredictions(int *predictions) {

    printf("\n\n\n------------- PREDICTIONS -------------\n\n");
    int currentMax = 1000;

    for (int p = 0; p < 10; p++) {
        int highest = 0;
        int value = 0;

        for (int i = 0; i < MAX_CLASS - 1; i++) {
            if (predictions[i] > value && predictions[i] < currentMax) {
                highest = i;
                value = predictions[i];
            }
        }
        std::cout << "[" << highest << "]" << labels[highest] << ": " << value << " ("
                  << (float) ((float) value / (float) TREE_NUM) * 100.0 << "%)" << std::endl;
        currentMax = value;
    }
    printf("\n");


}

/**
 *
 * Sets the initial values for device variables
 *
 * @param currentNode
 * @param maxNode
 * @param furthestPossibleValidNode
 * @param sharedY
 * @param d_y
 * @param sharedMask
 * @param d_mask
 * @param forest
 * @param nodeMapping
 * @param tid
 */
__device__ void
initSharedValues(int *currentNode, int *maxNode, int *furthestPossibleValidNode,
                 int *sharedY, const int *d_y, int *sharedMask, int *d_mask, int *forest, int *nodeMapping, int tid, bool crossValidation) {

    if (threadIdx.x == 0 && threadIdx.y == 0) {
        *currentNode = 0;
        *maxNode = (1 << (MAX_DEPTH + 1)) - 1;
        *furthestPossibleValidNode = 4000;
    }

    // Put y into shared memory
    for (int i = tid; i < SAMPLES; i += blockDim.x * blockDim.y) {
        sharedY[i] = d_y[i];
        if(crossValidation) {
            sharedMask[i] = d_mask[i];
        } else {
            sharedMask[i] = i;
        }
        nodeMapping[i] = 0;
    }

    __syncthreads();
}


/**
 *
 * Resetting values to a default state between each node in each tree
 *
 * @param sampleSize
 * @param bestFeature
 * @param isLeaf
 * @param classDiscovered
 * @param classCount
 * @param branchClassCounts
 * @param tid
 */
__device__ void
reinitSharedValues(int *sampleSize, int *bestFeature, bool *isLeaf, int *classDiscovered, int *classCount,
                   int *branchClassCounts, int tid) {
    if (threadIdx.x == 0 && threadIdx.y == 0) {
        *sampleSize = 0;
        *bestFeature = 0;
        *isLeaf = false;
        *classDiscovered = 0;
    }

    // Reset global class counts
    if (tid < MAX_CLASS) {
        classCount[tid] = 0;
    }

    // Reset local class counts
    for (int i = threadIdx.x; i < MAX_CLASS * FEATURE_SQRT; i += blockDim.x) {
        branchClassCounts[i] = 0;
    }

    __syncthreads();
}


/**
 *
 * Resetting branches to default values for each node
 *
 * @param sharedLeftBranchEntropy
 * @param sharedRightBranchEntropy
 * @param sharedBranchSums
 * @param parentEntropy
 */
__device__ void resetBranches(float *sharedLeftBranchEntropy, float *sharedRightBranchEntropy, int *sharedBranchSums,
                              float *parentEntropy) {
    if (threadIdx.x == 0) {
        sharedLeftBranchEntropy[threadIdx.y] = 0;
        sharedRightBranchEntropy[threadIdx.y] = 0;
        sharedBranchSums[threadIdx.y] = 0;
        *parentEntropy = 0;
    }
    __syncthreads();
}


/**
 *
 * Go to the next node
 *
 * @param currentNode
 */
__device__ void nextNode(int *currentNode) {
    if (threadIdx.x == 0 && threadIdx.y == 0) {
        atomicAdd(currentNode, 1);
    }
    __syncthreads();
}


/**
 *
 * Get a random feature / column from the dataset to test for best split
 *
 * @param currentNode
 * @param state
 * @return
 */

__device__ int getRandomFeature(int currentNode, curandState *state) {
    return curand_uniform(state) * (FEATURE_SIZE - 1);
}


/**
 *
 * A method to check if the current node is a leaf node
 *
 * @param isLeaf
 * @param classCount
 * @param classDiscovered
 * @param tid
 */
__device__ void isCurrentNodeLeaf(bool *isLeaf, int *classCount, int *classDiscovered, int tid) {

    // Check if the current node is a leaf node
    for (int i = tid; i < MAX_CLASS; i += blockDim.x * blockDim.y) {
        if (classCount[tid] > 0) {
            atomicAdd(classDiscovered, 1);
        }
    }

    if (threadIdx.x == 0 && threadIdx.y == 0 && *classDiscovered == 1) {
        *isLeaf = true;
    }

    __syncthreads();

}

/**
 *
 * Creates a leaf node
 *
 * @param forest
 * @param classCount
 * @param maxNode
 * @param currentNode
 * @param leafCount
 */
__device__ void createLeafNode(int *forest, int *classCount, int maxNode, int currentNode, int *leafCount) {
    if (threadIdx.x == 0 && threadIdx.y == 0) {

        // Disable children if they exist
        if (currentNode * 2 + 2 <= maxNode) {
            forest[maxNode * blockIdx.x + (currentNode * 2 + 1)] = -1;
            forest[maxNode * blockIdx.x + (currentNode * 2 + 2)] = -1;
        }

        // Set leaf value
        int highestCount = 0;
        int highestClass = -1;
        for (int i = 0; i < MAX_CLASS; i++) {

            if (classCount[i] > highestCount) {
                highestCount = classCount[i];
                highestClass = i;
            }

        }

        if (highestClass >= 0) {
            forest[maxNode * blockIdx.x + currentNode] = -highestClass;
        } else {
            forest[maxNode * blockIdx.x + currentNode] = -69;
        }

        *leafCount = *leafCount + 1;
    }
    __syncthreads();
}


/**
 *
 * Exports the class count to global memory so it can be used to calculate
 * prediction probabilities later
 *
 * @param classCount
 * @param leafSamples
 * @param tid
 */
__device__ void setLeafClassCount(int *classCount, int *leafSamples, int tid) {

    if (tid < MAX_CLASS) {
        leafSamples[tid] = classCount[tid];
    }

    __syncthreads();
}


/**
 *
 * Calculates entropy for a given branch
 *
 * @param sampleSize
 * @param branchSize
 * @param tid
 * @param classCount
 * @param branchClassCounts
 * @param sharedLeftBranchEntropy
 * @param sharedRightBranchEntropy
 * @param parentEntropy
 */
__device__ void calculateEntropy(int sampleSize, int branchSize, int tid, int *classCount, int *branchClassCounts,
                                 float *sharedLeftBranchEntropy, float *sharedRightBranchEntropy,
                                 float *parentEntropy) {

    // Calculate parent entropy

    if (tid < MAX_CLASS) {
        float div = (float) classCount[tid] / (float) sampleSize;
        if (div > 0) {
            atomicAdd(parentEntropy, -(div * __log2f(div)));
        }
    }


    float leftSum = 0;
    float rightSum = 0;
    int rightBranchSize = sampleSize - branchSize;

    // Calculate partial left- and right-entropy
    if (branchSize > 0 && rightBranchSize > 0) {
        for (int i = threadIdx.x; i < MAX_CLASS; i += blockDim.x) {

            float localClassCount = (float) (branchClassCounts + MAX_CLASS * threadIdx.y)[i];
            float divLeft = localClassCount / (float) branchSize;
            float divRight =
                    (float) (classCount[i] - localClassCount) /
                    (float) rightBranchSize;


            if (divLeft > 0) {
                leftSum += divLeft * __log2f(divLeft);
            }

            if (divRight > 0) {
                rightSum += divRight * __log2f(divRight);
            }
        }

    }
    __syncthreads();

    // Add the partial entropy respectively
    atomicAdd(&sharedLeftBranchEntropy[threadIdx.y], -leftSum);
    atomicAdd(&sharedRightBranchEntropy[threadIdx.y], -rightSum);

    __syncthreads();
}


/**
 *
 * Calculates the weighted entropy of the left and right branches
 *
 * @param parentEntropy
 * @param sampleSize
 * @param branchSize
 * @param sharedLeftBranchEntropy
 * @param sharedRightBranchEntropy
 */
__device__ void
calculateInformationGain(float parentEntropy, int sampleSize, int branchSize, float *sharedLeftBranchEntropy,
                         float *sharedRightBranchEntropy) {
    // Calculate information gain for current feature
    if (threadIdx.x == 0) {
        int rightBranchSize = sampleSize - branchSize;

        if (rightBranchSize > 0 && branchSize > 0) {

            sharedLeftBranchEntropy[threadIdx.y] = parentEntropy -
                                                   (__fdividef((float) branchSize, (float) sampleSize) *
                                                    sharedLeftBranchEntropy[threadIdx.y] +
                                                    __fdividef((float) rightBranchSize, (float) sampleSize) *
                                                    sharedRightBranchEntropy[threadIdx.y]);
        } else {
            sharedLeftBranchEntropy[threadIdx.y] = -1000;
        }
    }

    __syncthreads();
}


/**
 *
 * Finds the feature with highest information gain
 *
 * @param bestGain
 * @param bestFeature
 * @param informationGains
 */
__device__ void selectBestFeature(float *bestGain, int *bestFeature, float *informationGains) {

    if (threadIdx.x == 0 && threadIdx.y == 0) {

        *bestGain = informationGains[0];
        *bestFeature = 0;
        for (int i = 1; i < FEATURE_SQRT; i++) {
            if (informationGains[i] > *bestGain) {
                *bestGain = informationGains[i];
                *bestFeature = i;
            }
        }
    }

    __syncthreads();

}


/**
 *
 * Runs One vs. Rest logistic regression model prediction
 *
 * @param X
 * @param weights
 * @param bias
 */
__global__ void predictCalibrationModel(float *X, float *weights, float *bias) {
    int tid = threadIdx.x * blockDim.y + threadIdx.y;

    if(tid < MAX_CLASS && X[tid] > 0) {
        double prediction = (double) 1.0 / ((double)1.0 + (double)exp(-(X[tid] * weights[tid] + bias[tid])));

        if(tid == 9) {
            printf("%d[%f] = %f (%f, %f) -- more : %f and %f\n", tid, X[tid], prediction, weights[tid], bias[tid], ((double )1.0 + (double)exp(-(X[tid] * weights[tid] + bias[tid]))), (double)-(X[tid] * weights[tid] + bias[tid]));
        }

    }
}

__device__ float convergenceTreshold = 0.00001;

/**
 *
 * Fits the calibration layer in a One vs. Rest method
 *
 * @param X
 * @param y
 * @param mask
 * @param weights
 * @param globalBias
 */
__global__ void fitCalibrationModel(float *X, int *y, int* mask, float *weights, float *globalBias) {

    __shared__ int sharedY[SAMPLES];
    __shared__ float diff[SAMPLES];
    __shared__ float cost[SAMPLES];
    __shared__ float sharedWeight;
    __shared__ float bias;
    __shared__ float learningRate;
    __shared__ int currentEpoch;
    __shared__ int epochs;
    __shared__ float prevCost;
    __shared__ bool converge;

    int sampleSize = VALIDATION_FOLD_SIZE;

    int tid = threadIdx.x * blockDim.y + threadIdx.y;


    if(tid < sampleSize) {
        int index = tid;
        sharedY[index] = y[mask[blockIdx.y * VALIDATION_FOLD_SIZE + index]] == blockIdx.x;
    }


    if(threadIdx.x == 0 && threadIdx.y == 0) {
        epochs = 10000;
        currentEpoch = 0;
        learningRate = 0.1;
        bias = 1;
        sharedWeight = 0.1;
        converge = false;
    }

    while(currentEpoch < epochs) {
        if(converge) {
            break;
        }

     __shared__   if(tid < sampleSize) {
            if(blockIdx.x == 9) {
                //printf("[%f, %d]\n", X[1 * VALIDATION_FOLD_SIZE+tid+blockIdx.x], sharedY[tid]);
            }
            float prediction = 1.0 / (1.0 + exp(-(X[blockIdx.y * VALIDATION_FOLD_SIZE+tid+blockIdx.x] * sharedWeight + bias)));

            float epsilon = 1e-7;
            prediction = max(min(prediction, 1.0f - epsilon), epsilon);

            diff[tid] = prediction - sharedY[tid];

            cost[tid] = -(sharedY[tid] * log(prediction) + (1 - sharedY[tid]) * log(1 - prediction));

            if(currentEpoch == 0 && blockIdx.x == 9) {
                //printf("%f,", X[tid]);
            }

        }

        if(tid == 0) {

            float sum = 0.0;
            float diffSum = 0.0;
            float costSum = 0.0;

            for(int i = 0; i < sampleSize; i++) {
                sum += diff[i] * X[blockIdx.y * VALIDATION_FOLD_SIZE+i+blockIdx.x];
                diffSum += diff[i];
                costSum += cost[i];
            }

            sharedWeight -= learningRate * ((1.0 / (float) sampleSize) * sum);
            bias -= learningRate * ((1.0 / (float) sampleSize) * diffSum);

            if(currentEpoch > 0) {
                if(fabs(prevCost - costSum) < convergenceTreshold) {
                    converge = true;
                }
            }

            prevCost = costSum;
            currentEpoch++;
        }

        __syncthreads();
    }

    if(tid == 0) {
        if(blockIdx.x == 9) {
            printf("Stopped at iteration %d with loss=%f\n", currentEpoch, prevCost);
        }
        weights[blockIdx.x] = sharedWeight;
        globalBias[blockIdx.x] = bias;
    }



}


/**
 *
 * Runs prediction on a decision tree
 *
 * @param forest
 * @param d_x
 * @param predictions
 * @param leafSamples
 * @param sharedPredictionProba
 * @param treeSize
 * @param tree
 * @param validation
 */
__device__ void predictTree(const int* forest, const int* d_x, int* predictions, const int* leafSamples, float* sharedPredictionProba, int treeSize, int tree, bool validation) {
    int currentNode = 0;
    int currentValue = 0;

    while (currentNode < treeSize) {
        currentValue = forest[tree * treeSize + currentNode];
        if (currentValue >= 0) {
            if (d_x[currentValue] == 1) {
                currentNode = currentNode * 2 + 2;
            } else {
                currentNode = currentNode * 2 + 1;
            }
        } else {
            int leafValue = -currentValue;

            if (leafValue < MAX_CLASS) {

                if(!validation) {
                    atomicAdd(&predictions[leafValue], 1);
                }

                int samples = 0;
                for (int i = 0; i < MAX_CLASS; i++) {
                    samples += leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i];
                }

                for (int i = 0; i < MAX_CLASS; i++) {
                    if (leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i] >= 0 && samples > 0) {
                        float proba =
                                (float) leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i] / (float) samples;
                        atomicAdd(&sharedPredictionProba[leafValue], proba);
                    }
                }
            }
            break;
        }
    }
}


/**
 *
 * Validation prediction probabilities for a decision tree
 *
 * @param forest
 * @param d_x
 * @param predictions
 * @param leafSamples
 * @param validationProba
 * @param treeSize
 * @param tree
 * @param validation
 */
__device__ void validateTree(const int* forest, const int* d_x, int* predictions, const int* leafSamples, float* validationProba, int treeSize, int tree, bool validation) {
    int currentNode = 0;
    int currentValue = 0;

    while (currentNode < treeSize) {
        currentValue = forest[tree * treeSize + currentNode];
        if (currentValue > 0) {
            if (d_x[currentValue] == 1) {
                currentNode = currentNode * 2 + 2;
            } else {
                currentNode = currentNode * 2 + 1;
            }
        } else {
            int leafValue = -currentValue;


            if (leafValue >= 0) {

                int samples = 0;
                for (int i = 0; i < MAX_CLASS; i++) {
                    samples += leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i];
                }


                if(samples > 0) {
                    for (int i = 0; i < MAX_CLASS; i++) {
                        if (leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i] > 0) {
                            float proba =
                                    (float) leafSamples[tree * LEAF_NODES_PER_TREE + MAX_CLASS * currentNode + i] / (float) samples;
                            atomicAdd(&validationProba[i], proba);
                        }
                    }
                }



            }
            break;
        }
    }
}

/**
 *
 * Run prediction on every decision tree, producing hard votes and prediction probabilities.
 * Also produces validation probabilities if it is in a cross validation fold.
 *
 * @param forest
 * @param leafSamples
 * @param d_x
 * @param sample
 * @param validationFolds
 * @param treeSize
 * @param predictions
 * @param predictionProba
 * @param validationProba
 * @param validationFoldIndex
 * @param fold
 * @param probaTransfer
 */
__global__ void predictForest(const int *forest, const int *leafSamples, const int *d_x, const int *sample, int* validationFolds, int treeSize, int *predictions,
                              float *predictionProba, float *validationProba, int validationFoldIndex, int fold, float* probaTransfer) {

    __shared__ float sharedPredictionProba[MAX_CLASS];
    __shared__ int count;

    int tid = threadIdx.x;
    int tree = tid;
    int validationSample = 0;


    if(tid == 0) {
        count = 0;
    }

    if (tid < MAX_CLASS) {
        sharedPredictionProba[tid] = 0;
    }

    __syncthreads();

    while (tree < TREE_NUM) {
        validationSample = 0;
        predictTree(forest, sample, predictions, leafSamples, sharedPredictionProba, treeSize, tree, false);

        while(validationSample < VALIDATION_FOLD_SIZE) {
            validateTree(forest, d_x + (validationFolds[validationSample]) * FEATURE_SIZE, predictions, leafSamples, validationProba + (fold*VALIDATION_FOLD_SIZE + validationSample) * MAX_CLASS, treeSize, tree, true);
            validationSample++;
        }
        tree += blockDim.x * blockDim.y;
        atomicAdd(&count, 1);
    }

    __syncthreads();

    if (tid < MAX_CLASS) {
        predictionProba[tid] = sharedPredictionProba[tid] / (float) TREE_NUM;
        validationSample = 0;
        while(validationSample < VALIDATION_FOLD_SIZE) {
            (validationProba + (fold*VALIDATION_FOLD_SIZE + validationSample) * MAX_CLASS)[tid] /= (float) count;
            validationSample++;
        }

    }

    __syncthreads();
    if(tid == 0) {
        validationSample = 0;
        while (validationSample < VALIDATION_FOLD_SIZE) {


            int highestIndex = -1;
            float highestProba = 0;

            for (int i = 0; i < MAX_CLASS; i++) {

                float value = (validationProba + (fold * VALIDATION_FOLD_SIZE + validationSample) * MAX_CLASS)[i];
                if (value > highestProba) {
                    highestProba = value;
                    highestIndex = i;
                }
            }

            probaTransfer[validationSample*2] = highestIndex;
            probaTransfer[validationSample*2+1] = highestProba;
            validationSample++;
        }
    }


}



/**
 *
 * CUDA kernel for building the random forest
 *
 * @param X
 * @param y
 * @param mask
 * @param branchClassCounts
 * @param leafSamples
 * @param forest
 * @param samples
 * @param featureSize
 * @param seed
 * @param validationFoldIndex
 * @param crossValidate
 */
__global__ void
buildForest(int *X, int *y, int *mask, int *branchClassCounts, int *leafSamples, int *forest, int samples,
            int featureSize, int seed, int validationFoldIndex, bool crossValidate) {

    extern __shared__ int sharedData[];
    int* sharedMask = sharedData;
    int* sharedY = sharedData + SAMPLES;
    int* nodeMapping = sharedData + SAMPLES * 2;
    int* classCount = sharedData + SAMPLES * 3;
    int* features = sharedData + SAMPLES * 3 + MAX_CLASS;
    int* sharedBranchSums = sharedData + SAMPLES * 3 + MAX_CLASS + FEATURE_SQRT;
    float* sharedRightBranchEntropy = (float *) sharedData + SAMPLES * 3 + MAX_CLASS + FEATURE_SQRT*2;
    float* sharedLeftBranchEntropy = (float *) sharedData + SAMPLES * 3 + MAX_CLASS + FEATURE_SQRT*3;

    __shared__ int currentNode;
    __shared__ int sampleSize;
    __shared__ int classDiscovered;
    __shared__ float parentEntropy;
    __shared__ int bestFeature;
    __shared__ int maxNode;
    __shared__ int furthestPossibleValidNode;
    __shared__ float bestGain;
    __shared__ bool isLeaf;
    __shared__ int leafCount;


    if(blockIdx.x == 0 && blockIdx.y == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        FEATURE_SIZE = featureSize;
    }

    int tid = threadIdx.y * blockDim.x + threadIdx.x;

    __syncthreads();

    // Initialize values once
    initSharedValues(&currentNode, &maxNode, &furthestPossibleValidNode, sharedY, y, sharedMask, mask, forest,
                     nodeMapping, tid, crossValidate);

    if (threadIdx.x == 0 && threadIdx.y == 0) {
        leafCount = 0;
    }

    curandState state;
    curand_init(seed, blockIdx.x + blockIdx.y * gridDim.x * (blockDim.x * blockDim.y) + threadIdx.x +
                      blockDim.x * threadIdx.y, 0, &state);

    __syncthreads();

    while (currentNode < maxNode) {

        reinitSharedValues(&sampleSize, &bestFeature, &isLeaf, &classDiscovered, classCount,
                           branchClassCounts + MAX_CLASS * FEATURE_SQRT * blockIdx.x, tid);

        if (forest[maxNode * blockIdx.x + currentNode] != -1) {
            resetBranches(sharedLeftBranchEntropy, sharedRightBranchEntropy, sharedBranchSums, &parentEntropy);

            // Grab feature
            features[threadIdx.y] = getRandomFeature(currentNode, &state);
            __syncthreads();

            // Count occurrences of the different classes for parent and each feature
            int i = threadIdx.x;
            int branchSum = 0;

            while (i < SAMPLES) {
                int maskIndex = sharedMask[i];
                if (nodeMapping[maskIndex] == currentNode) {
                    if (threadIdx.y == 0) {
                        atomicAdd(&classCount[sharedY[maskIndex]], 1);
                        atomicAdd(&sampleSize, 1);
                    }

                    if (X[maskIndex * FEATURE_SIZE + features[threadIdx.y]] == 0) {
                        atomicAdd(&(branchClassCounts + blockIdx.x * MAX_CLASS * FEATURE_SQRT +
                                    MAX_CLASS * threadIdx.y)[sharedY[maskIndex]], 1);
                        branchSum++;
                    }
                }

                if (i >= validationFoldIndex && i < validationFoldIndex + (SAMPLES / 5)) {
                    i += SAMPLES / 5;
                } else {
                    i += blockDim.x;
                }

            }
            __syncthreads();

            // Calculate the branch size
            atomicAdd(&sharedBranchSums[threadIdx.y], branchSum);

            __syncthreads();

            // Check if this node should be a leaf
            isCurrentNodeLeaf(&isLeaf, classCount, &classDiscovered, tid);

            int leftChild = currentNode * 2 + 1;


            // Stopping criteria
            if (isLeaf || sampleSize <= 0 || sampleSize < MIN_SAMPLES_SPLIT || leftChild + 1 >= maxNode - 1 ||
                sampleSize <= MIN_SAMPLES_PER_LEAF) {

                // Create leaf node and store the class samples
                createLeafNode(forest, classCount, maxNode, currentNode, &leafCount);
                setLeafClassCount(classCount, leafSamples + blockIdx.x * LEAF_NODES_PER_TREE + currentNode * MAX_CLASS, tid);
            } else {

                // Select the feature with the highest information gain
                calculateEntropy(sampleSize, sharedBranchSums[threadIdx.y], tid, classCount,
                                 branchClassCounts + MAX_CLASS * FEATURE_SQRT * blockIdx.x,
                                 sharedLeftBranchEntropy, sharedRightBranchEntropy, &parentEntropy);
                calculateInformationGain(parentEntropy, sampleSize, sharedBranchSums[threadIdx.y],
                                         sharedLeftBranchEntropy, sharedRightBranchEntropy);

                selectBestFeature(&bestGain, &bestFeature, sharedLeftBranchEntropy);


                // If positive information gain, split the tree
                if (bestGain > IG_THRESHOlD) {

                    // Assign the feature value to the current node
                    if (threadIdx.x == 0 && threadIdx.y == 0) {
                        forest[maxNode * blockIdx.x + currentNode] = features[bestFeature];
                    }
                    __syncthreads();

                    // Update node mapping
                    int index = tid;
                    while (index < SAMPLES) {

                        int maskIndex = sharedMask[index];
                        if (nodeMapping[maskIndex] == currentNode) {
                            nodeMapping[maskIndex] = leftChild;
                            if (X[maskIndex * FEATURE_SIZE + features[bestFeature]] == 1) {
                                nodeMapping[maskIndex]++;
                            }
                        }

                        if (index >= validationFoldIndex && index < validationFoldIndex + SAMPLES / 5) {
                            index += SAMPLES / 5;
                        } else {
                            index += blockDim.x * blockDim.y;
                        }
                    }

                } else {
                    // If no positive gain can be found, create a leaf node and prevent further splitting
                    createLeafNode(forest, classCount, maxNode, currentNode, &leafCount);
                    setLeafClassCount(classCount, leafSamples + blockIdx.x * LEAF_NODES_PER_TREE + currentNode * MAX_CLASS, tid);
                }
            }

        }
        __syncthreads();

        nextNode(&currentNode);
    }
}


/**
 *
 * Allocate memory on GPU, copy data from host and set memory areas
 *
 * @param d_x
 * @param h_x
 * @param d_y
 * @param h_y
 * @param d_mask
 * @param h_mask
 * @param d_trainingFolds
 * @param h_trainingFolds
 * @param d_validationFolds
 * @param h_validationFolds
 * @param d_leafSamples
 * @param d_branchClassCounts
 * @param d_forest
 * @param d_predictions
 * @param d_predictionProba
 * @param d_validationProba
 * @param d_x_test
 * @param d_weights
 * @param d_bias
 * @param h_test
 * @param sampleSize
 * @param featureSize
 * @param featuresUsed
 * @param maxLeafNodes
 * @param treeSize
 */
void allocateDeviceMemory(int **d_x, int *h_x, int **d_y, int *h_y, int** d_mask, int* h_mask, int **d_trainingFolds, int* h_trainingFolds, int** d_validationFolds, int* h_validationFolds, int **d_leafSamples,
                          int **d_branchClassCounts, int **d_forest, int **d_predictions, float **d_predictionProba,
                          float **d_validationProba,
                          int **d_x_test, float **d_weights, float **d_bias, int *h_test, int sampleSize, int featureSize, int featuresUsed,
                          int maxLeafNodes, int treeSize) {

     printf("Feature size in allocateDeviceMemory : %d\n", featureSize);
    // Allocate memory on the GPU device
    cudaMalloc(d_x, sizeof(int) * SAMPLES * featureSize);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_x)");

    cudaMalloc(d_y, sizeof(int) * SAMPLES);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_y)");

    cudaMalloc(d_mask, sizeof(int) * SAMPLES *5);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_mask)");

    cudaMalloc(d_trainingFolds, sizeof(int) * TRAINING_FOLD_SIZE * 5);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_trainingFolds)");

    cudaMalloc(d_validationFolds, sizeof(int) * VALIDATION_FOLD_SIZE * 5);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_validationFolds)");

    cudaMalloc(d_leafSamples, sizeof(int) * MAX_CLASS * maxLeafNodes * TREE_NUM);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_leafSamples)");

    cudaMalloc(d_branchClassCounts, sizeof(int) * MAX_CLASS * featuresUsed * TREE_NUM);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_branchClassCounts)");

    cudaMalloc(d_forest, sizeof(int) * treeSize * TREE_NUM);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_forest)");

    cudaMalloc(d_predictions, sizeof(int) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_predictions)");

    cudaMalloc(d_predictionProba, sizeof(float) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_predictionProba)");

    cudaMalloc(d_validationProba, sizeof(float) * MAX_CLASS * SAMPLES);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_validationProba)");

    cudaMalloc(d_x_test, sizeof(int) * featureSize);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_x_test)");

    cudaMalloc(d_weights, sizeof(float) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_weights)");

    cudaMalloc(d_bias, sizeof(float) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory (d_bias)");

    // Copy data from host to device
    cudaMemcpy(*d_x, h_x, sizeof(int) * SAMPLES * featureSize, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_x)");

    cudaMemcpy(*d_y, h_y, sizeof(int) * SAMPLES, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_y)");

    cudaMemcpy(*d_mask, h_mask, sizeof(int) * SAMPLES *5, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_mask)");

    cudaMemcpy(*d_trainingFolds, h_trainingFolds, sizeof(int) * TRAINING_FOLD_SIZE * 5, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_trainingFolds)");

    cudaMemcpy(*d_validationFolds, h_validationFolds, sizeof(int) * VALIDATION_FOLD_SIZE * 5, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_validationFolds)");

    cudaMemcpy(*d_x_test, h_test, sizeof(int) * featureSize, cudaMemcpyHostToDevice);
    checkCudaError(cudaStatus, "allocateDeviceMemory[COPY] (d_x_test)");

    // Initialize the calibration weights
    cudaMemset(*d_weights, 0.01, sizeof(float) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory[MEMSET] (d_weights)");
    cudaMemset(*d_bias, 0, sizeof(float) * MAX_CLASS);
    checkCudaError(cudaStatus, "allocateDeviceMemory[MEMSET] (d_bias)");
    cudaMemset(*d_validationProba, 0, sizeof(float) * MAX_CLASS * SAMPLES);
    checkCudaError(cudaStatus, "allocateDeviceMemory[MEMSET] (d_validationProba)");

}

/**
 *
 * Free memory on the GPU device
 *
 * @param d_x
 * @param d_y
 * @param d_mask
 * @param d_leafSamples
 * @param d_branchClassCounts
 * @param d_tree
 */
void freeDeviceMemory(int *d_x, int *d_y, int *d_mask, int *d_leafSamples, int *d_branchClassCounts, int *d_tree) {
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_mask);
    cudaFree(d_leafSamples);
    cudaFree(d_branchClassCounts);
    cudaFree(d_tree);
}


/**
 *
 * Fits a random forest on the entire dataset
 *
 * @param d_tree
 * @param d_leafSamples
 * @param d_predictions
 * @param d_predictionProba
 * @param d_validationProba
 * @param d_branchClassCounts
 * @param d_mask
 * @param d_x
 * @param d_y
 * @param d_sample
 * @param h_predictions
 * @param seed
 * @param treeSize
 * @param maxLeafNodes
 * @param memory
 * @param featureSize
 */
void fitBaseRF(int **d_tree, int **d_leafSamples, int **d_predictions, float **d_predictionProba, float **d_validationProba, int **d_branchClassCounts, int **d_mask, int **d_x, int **d_y, int **d_sample, int **h_predictions, int seed, int treeSize, int maxLeafNodes, int memory, int featureSize) {
    // Reset values
    cudaMemset(*d_tree, 0, sizeof(int) * treeSize * TREE_NUM);
    cudaMemset(*d_leafSamples, 0, sizeof(int) * MAX_CLASS * maxLeafNodes * TREE_NUM);
    cudaMemset(*d_predictions, 0, sizeof(int) * MAX_CLASS);
    cudaMemset(*d_predictionProba, 0, sizeof(float) * MAX_CLASS);
    cudaMemset(*d_validationProba, 0, sizeof(float) * MAX_CLASS * SAMPLES);

    printf("Training Random Forest...\n");

    // Build a Random Forest
    buildForest<<<gridSize, blockSize, memory>>>(*d_x, *d_y, *d_mask, *d_branchClassCounts, *d_leafSamples, *d_tree, SAMPLES,
                                         featureSize, seed, SAMPLES+1, false);
    cudaDeviceSynchronize();
    checkCudaError(cudaStatus, "buildForest");

    // Predict the test sample and the validation set
    predictForest<<<predictionGrid, predictionBlock>>>(*d_tree, *d_leafSamples, *d_x, *d_sample, *d_mask, treeSize, *d_predictions,
                                                       *d_predictionProba, *d_validationProba, SAMPLES+1, 0, *d_validationProba);
    cudaDeviceSynchronize();
    checkCudaError(cudaStatus, "predictForest");

    // Print the predictions
    cudaMemcpy(*h_predictions, *d_predictions, sizeof(int) * MAX_CLASS, cudaMemcpyDeviceToHost);
    printPredictions(*h_predictions);
}

struct Node{
    Node* left;
    Node* right;
    int feature;
    int value;
    bool isLeaf;
};

void freeBinaryForest(Node* node) {

    if(!node->isLeaf) {
        freeBinaryForest(node->left);
        freeBinaryForest(node->right);
    }

    free(node);

}

Node* createBinaryTree(int* forest, int forestSize, int index) {


    Node* node = (Node *) malloc(sizeof(Node));

    if(forest[index] >= 0) {
        node->feature = forest[index];

        node->left = createBinaryTree(forest, forestSize, index*2+1);
        node->right = createBinaryTree(forest, forestSize, index*2+2);
    } else {
        node->value = -forest[index];
        node->isLeaf = true;
    }


    return node;
}

void buildVisualForest(int* forest, int treeSize) {

    //Node* root = createBinaryTree(forest, treeSize, 0);


}

/**
 *
 * The fit method of the classifier.
 * Performs Cross validation and fitting of the calibration layer.
 * Then trains a random forest on the entire dataset and gets it's predictions.
 * The predictions are adjusted by the calibration layer.
 *
 * @param X
 * @param y
 * @param sample
 * @param foldBitMask
 * @param trainingFolds
 * @param validationFolds
 * @param sampleSize
 * @param featureSize
 * @param seed
 * @return
 */
ClassifierResult CudaCalibratedClassifierCV::fit(int *X, int *y, int *sample, int* foldBitMask, int* trainingFolds, int* validationFolds, int sampleSize, int featureSize, int seed) {

    // TODO : Check if the prediction percentages are correct (Probably solved now)
    // TODO : Error checks
    // TODO : Compare calibration results
    // TODO : Implement Bootstrapping
    // TODO : Implement Out-of-bag error
    // TODO : Test edge cases (MAX_CLASS, TREE_NUM = 1, FEATURES, ...)

    // Host pointers
    int nFolds = 5;
    int treeSize = (1 << (MAX_DEPTH + 1)) - 1;
    int *h_predictions = (int *) malloc(sizeof(int) * MAX_CLASS);

    // Be wary : used to be Max depth + 1
    int maxLeafNodes = (1 << (MAX_DEPTH));
    int validationFoldIndex = 0;
    int foldSize = SAMPLES / nFolds;
    this->featureSize = featureSize;


    // Device pointers
    int *d_x;
    int *d_y;
    int *d_mask;
    int *d_trainingFolds;
    int *d_validationFolds;
    int *d_leafSamples;
    int *d_branchClassCounts;
    int *d_tree;
    int *d_predictions;
    float *d_predictionProba;
    float *d_validationProba;
    float *d_weights;
    float *d_bias;
    int *d_sample;
    float *d_probaTransfer;
    cudaMalloc(&d_probaTransfer, sizeof(float) * VALIDATION_FOLD_SIZE *2);


    int totalSharedMemory = SAMPLES * sizeof(int) * 3 + MAX_CLASS * sizeof(int)
                            + FEATURE_SQRT * sizeof(int) * 4;
    printf("Running with dimension : (%d, %d) and total shared memory : %d\n", sampleSize, featureSize, totalSharedMemory);


    allocateDeviceMemory(&d_x, X, &d_y, y, &d_mask, foldBitMask,&d_trainingFolds, trainingFolds, &d_validationFolds, validationFolds, &d_leafSamples, &d_branchClassCounts, &d_tree,
                         &d_predictions, &d_predictionProba, &d_validationProba, &d_sample, &d_weights, &d_bias, sample, sampleSize, featureSize, FEATURE_SQRT,
                         maxLeafNodes, treeSize);
    checkCudaError(cudaStatus, "allocateDeviceMemory");

    int tempx = 0;
    printf("\n\n---------- Running 5-fold cross-validation ---------- \n\n");
    for (int i = 0; i < nFolds; i++) {
        validationFoldIndex = SAMPLES - foldSize;
        // Reset values
        cudaMemset(d_tree, 0, sizeof(int) * treeSize * TREE_NUM);
        cudaMemset(d_leafSamples, 0, sizeof(int) * MAX_CLASS * maxLeafNodes * TREE_NUM);
        cudaMemset(d_predictions, 0, sizeof(int) * MAX_CLASS);
        cudaMemset(d_predictionProba, 0, sizeof(float) * MAX_CLASS);


        printf("Running fold %d\n", i+1);

        // Build a Random Forest
        buildForest<<<gridSize, blockSize, totalSharedMemory>>>(d_x, d_y, d_trainingFolds + TRAINING_FOLD_SIZE * i, d_branchClassCounts, d_leafSamples, d_tree, sampleSize,
                                             featureSize, seed, validationFoldIndex, true);
        cudaDeviceSynchronize();
        checkCudaError(cudaStatus, "buildForest");

        // Predict the test sample and the validation set
        predictForest<<<predictionGrid, predictionBlock>>>(d_tree, d_leafSamples, d_x, d_sample, d_validationFolds + VALIDATION_FOLD_SIZE*i, treeSize, d_predictions,
                                                           d_predictionProba, d_validationProba, validationFoldIndex, tempx, d_probaTransfer);
        cudaDeviceSynchronize();
        checkCudaError(cudaStatus, "predictForest");

        // Print the predictions
        //cudaMemcpy(h_predictions, d_predictions, sizeof(int) * MAX_CLASS, cudaMemcpyDeviceToHost);
        //printPredictions(h_predictions);
        tempx++;

    }

    float *temp = (float *) malloc(sizeof(float) * MAX_CLASS * SAMPLES);
    cudaMemcpy(temp, d_validationProba, sizeof(float) * MAX_CLASS * SAMPLES, cudaMemcpyDeviceToHost);

    float *yup = (float *) malloc(sizeof(float) * VALIDATION_FOLD_SIZE * 2);
    cudaMemcpy(yup, d_probaTransfer, sizeof(float) * VALIDATION_FOLD_SIZE * 2, cudaMemcpyDeviceToHost);

    printf("Validation probabilities : \n");

    for(int s = 0; s < 5; s++) {
        float total = 0.0;

        for(int i = 0; i < MAX_CLASS; i++) {
            printf("%f,", temp[s*MAX_CLASS + i]);
            total += temp[s*MAX_CLASS+i];
        }
        printf("Total : %f\n", total);
    }
    printf("\n----------\n");


    // Train the calibration models
    fitCalibrationModel<<<calibrationGrid, calibrationBlock>>>(d_validationProba, d_y, d_validationFolds, d_weights, d_bias);
    cudaDeviceSynchronize();
    checkCudaError(cudaStatus, "fitCalibrationModel");

    fitBaseRF(&d_tree, &d_leafSamples, &d_predictions, &d_predictionProba, &d_validationProba, &d_branchClassCounts, &d_mask, &d_x, &d_y, &d_sample, &h_predictions, seed, treeSize, maxLeafNodes, totalSharedMemory, featureSize);

    float *h_proba = (float *) malloc(sizeof(float) * MAX_CLASS);
    cudaMemcpy(h_proba, d_predictionProba, sizeof(float) * MAX_CLASS, cudaMemcpyDeviceToHost);

    ClassifierResult result = ClassifierResult();
    result.proba = h_proba;


    printf("\n\n------------------- Calibrated scores -------------------\n");

    dim3 predictCalBlock(MAX_CLASS, 1);
    predictCalibrationModel<<<predictionGrid, predictCalBlock>>>(d_predictionProba, d_weights, d_bias);
    cudaDeviceSynchronize();

    cudaMemcpy(h_proba, d_predictionProba, sizeof(float) * MAX_CLASS, cudaMemcpyDeviceToHost);

    printf("Prediction probabilities : \n");
    for(int i = 0; i < MAX_CLASS; i++) {
        printf("%f,", h_proba[i]);
    }
    printf("\n\n");


    //freeDeviceMemory(d_x, d_y, d_mask, d_leafSamples, d_branchClassCounts, d_tree);
    checkCudaError(cudaStatus, "freeDeviceMemory");
    printf("Freed device memory\n");

    result.predictions = h_predictions;
    result.calibratedProba = h_proba;

    return result;
}

void lrcFit(float *X, int *y, float *test, int sampleSize, int testSize, int featureSize){
    float* d_x;
    float* d_x_test;
    int* d_y;
    float* d_weights;
    float* d_bias;
    int* d_mask;

    int* h_mask = (int *) malloc(sizeof(int) * sampleSize);
    for(int i = 0; i < sampleSize; i++) {
        h_mask[i] = i;
    }



    cudaMalloc(&d_x, sizeof(float) * sampleSize * featureSize);
    cudaMemcpy(d_x, X, sizeof(float) * sampleSize * featureSize, cudaMemcpyHostToDevice);

    cudaMalloc(&d_y, sizeof(int) * sampleSize);
    cudaMemcpy(d_y, y, sizeof(int) * sampleSize, cudaMemcpyHostToDevice);

    cudaMalloc(&d_x_test, sizeof(float) * testSize * featureSize);
    cudaMemcpy(d_x_test, test, sizeof(float) * testSize * featureSize, cudaMemcpyHostToDevice);

    cudaMalloc(&d_mask, sizeof(int) * sampleSize);
    cudaMemcpy(d_mask, h_mask, sizeof(int) * sampleSize, cudaMemcpyHostToDevice);



    float* h_weights = (float *) malloc(sizeof(float) * MAX_CLASS);
    for(int i = 0; i < MAX_CLASS; i++) {
        h_weights[i] = (float) 0.1;
    }


    cudaMalloc(&d_weights, sizeof(float) * MAX_CLASS);
    cudaMalloc(&d_bias, sizeof(float) * MAX_CLASS);


    //cudaMemcpy(d_weights, h_weights, sizeof(float) * MAX_CLASS, cudaMemcpyHostToDevice);


    cudaMemset(d_bias, 0, sizeof(float) * MAX_CLASS);

    checkCudaError(cudaStatus, "cudaMemcpy");

    dim3 gridSize(1,1);
    dim3 blockSize(sampleSize,1);

    fitCalibrationModel<<<gridSize, blockSize>>>(d_x, d_y, d_mask, d_weights, d_bias);
    cudaDeviceSynchronize();

    dim3 blockSize2(1,1);
    predictCalibrationModel<<<gridSize, blockSize2>>>(d_x_test, d_weights, d_bias);
    cudaDeviceSynchronize();


    checkCudaError(cudaStatus, "cudaMemcpy");

    printf("Calibration model fit\n");
}

/**
 *
 * Constructor for the classifier class
 *
 * @param forestSize
 * @param maxDepth
 * @param maxClass
 * @param maxFeature
 * @param minSamplesLeaf
 */
CudaCalibratedClassifierCV::CudaCalibratedClassifierCV(int forestSize, int maxDepth, int maxClass, int maxFeature, int minSamplesLeaf){
    CudaCalibratedClassifierCV::forestSize = forestSize;
    CudaCalibratedClassifierCV::maxDepth = maxDepth;
    CudaCalibratedClassifierCV::maxClass = maxClass;
    CudaCalibratedClassifierCV::maxFeature = maxFeature;
    CudaCalibratedClassifierCV::minSampleLeaf = minSamplesLeaf;
}

/**
 * Prints properties for CUDA and classifier
 */
void CudaCalibratedClassifierCV::printProperties() {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);

    if (deviceCount == 0) {
        std::cerr << "No CUDA devices found!" << std::endl;
        return;
    }

    std::cout << "--------------- CUDA PROPERTIES ---------------" << std::endl;

    for (int device = 0; device < deviceCount; ++device) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);

        std::cout << "\tDevice " << device << ": " << deviceProp.name << std::endl;
        std::cout << "\tCompute capability: " << deviceProp.major << "." << deviceProp.minor << std::endl;
        std::cout << "\tTotal global memory: " << (deviceProp.totalGlobalMem / (1024 * 1024)) << " MB" << std::endl;
        std::cout << "\tMultiprocessors: " << deviceProp.multiProcessorCount << std::endl;
        std::cout << "\tClock rate: " << (deviceProp.clockRate / 1000) << " MHz" << std::endl;
        std::cout << "\tMemory clock rate: " << (deviceProp.memoryClockRate / 1000) << " MHz" << std::endl;
        std::cout << "\tMemory bus width: " << deviceProp.memoryBusWidth << " bits" << std::endl;
        // Add more properties as needed
    }
    std::cout << "-----------------------------------------------\n" << std::endl;

    std::cout << "--------------- CLASSIFIER PROPERTIES ---------------" << std::endl;
    std::cout << "\tFOREST SIZE : " << this->forestSize << std::endl;
    std::cout << "\tMAX DEPTH : " << this->maxDepth << std::endl;
    std::cout << "\tMAX FEATURE : " << this->maxFeature << std::endl;
    std::cout << "\tMIN SAMPLES PER LEAF : " << this->minSampleLeaf << std::endl;
    std::cout << "-----------------------------------------------------\n" << std::endl;
}


int CudaCalibratedClassifierCV::getFeatureSize() {
    return this->featureSize;
}

int CudaCalibratedClassifierCV::getForestSize() {
    return this->forestSize;
}

float CudaCalibratedClassifierCV::getOobError() {
    return this->oobError;
}