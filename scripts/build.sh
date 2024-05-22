rm lib/*
rm cudaclassifier*
rm -r build

/usr/local/cuda-12/bin/nvcc cuda/CudaCalibratedClassifierCV.cu  -Xcompiler -fPIC -lcudart -lcurand -lib -I common -o kernel -gencode arch=compute_75,code=sm_75

python cudaSetup.py build_ext -i