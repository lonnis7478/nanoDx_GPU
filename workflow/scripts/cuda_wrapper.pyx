cimport cython
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from cython.parallel cimport prange
import numpy as np

class classifier_result:
    def __init__(self, predictions, proba, calibrated_proba):
        self.predictions = predictions
        self.proba = proba
        self.calibrated_proba = calibrated_proba

    def predicted_class(self):
        return np.argmax(self.predictions)

cdef extern from "cuda/CudaCalibratedClassifierCV.h":

    cdef struct ClassifierResult:
            int* predictions
            float* proba
            float* calibratedProba

    cppclass CudaCalibratedClassifierCV:
        CudaCalibratedClassifierCV(int forestSize, int maxDepth, int maxClass, int maxFeature, int minSamplesLeaf)
        ClassifierResult fit(int* X, int* y, int* sample, int* foldBitmask,  int *trainingFolds, int* validationFolds, int samples, int featureSize, int seed)
        void printProperties()
        int getFeatureSize()
        int getForestSize()
        float getOobError()

cdef class CudaClassifier:
    cdef CudaCalibratedClassifierCV.CudaCalibratedClassifierCV *thisptr

    def __cinit__(self, n_ensemble=2000, max_depth=8, max_class=92, max_feature=96, min_samples_leaf=2):
        self.thisptr = new CudaCalibratedClassifierCV.CudaCalibratedClassifierCV(n_ensemble, max_depth, max_class, max_feature, min_samples_leaf)

    # Destructor
    def __dealloc__(self):
        del self.thisptr

    def print_cuda_properties(self):
        self.thisptr.printProperties()

    def get_feature_size(self):
        return self.thisptr.getFeatureSize()

    def get_forest_size(self):
        return self.thisptr.getForestSize()

    def get_oob_error(self):
        return self.thisptr.getOobError()

    def fit(self, int[::1] X, int[::1] y, int[::1] sample, int[::1] foldBitmask, int[::1] trainingFolds, int[::1] validationFolds, int sampleSize, int featureSize, int seed):
        cdef ClassifierResult result = self.thisptr.fit(&X[0], &y[0], &sample[0], &foldBitmask[0], &trainingFolds[0], &validationFolds[0], sampleSize, featureSize, seed)

        cr = classifier_result(np.array([result.predictions[i] for i in range(92)]),
                               np.array([result.proba[i] for i in range(92)]), np.array([result.calibratedProba[i] for i in range(92)]))

        return cr



