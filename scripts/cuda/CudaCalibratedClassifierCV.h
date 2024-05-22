#ifndef CUDACALIBRATEDCLASSIFIERCV_H
#define CUDACALIBRATEDCLASSIFIERCV_H


struct ClassifierResult {
    int *predictions;
    float *proba;
    float *calibratedProba;
};

class CudaCalibratedClassifierCV {

private:
    int forestSize;
    int featureSize;
    int maxDepth;
    int maxClass;
    int maxFeature;
    int minSampleLeaf;
    float oobError = 1.0;

public:

    // Constructor
    CudaCalibratedClassifierCV(int forestSize = 2000, int maxDepth = 8, int maxClass = 92, int maxFeature = 96,
                               int minSamplesLeaf = 1);

    ClassifierResult
    fit(int *X, int *y, int *sample, int *foldBitmask, int *trainingFolds, int *validationFolds, int samples,
        int featureSize, int seed);

    void printProperties();

    int getFeatureSize();

    float getOobError();

    int getForestSize();
};


void lrcFit(float *X, int *y, float *test, int sampleSize, int testSize, int featureSize);

#endif