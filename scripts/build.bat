del lib\\*.lib
del cudaext*
del cudaclassifier*
rmdir /s /q build

nvcc cuda/CudaCalibratedClassifierCV.cu -lcudart -lcurand -lib -I common -o lib/kernel.lib -gencode arch=compute_75,code=sm_75

python cudaSetup.py build_ext -i