#ifndef CYTHON_CUDA_SAMPLES_COMMON_H
#define CYTHON_CUDA_SAMPLES_COMMON_H
#include <cuda_runtime_api.h>
//#include <cuda.h>
#include <curand_kernel.h>
#include <string>
#include <iostream>


std::string labels[] = {"GBM- G34","DMG- K27","ATRT- SHH","CONTR- CEBM","GBM- MYCN","LGG- PA PF", "CONTR- REACT","LGG- PA MID","LGG- RGNT","MB- WNT",
                        "ATRT- MYC", "ATRT- TYR","LGG- PA/GG ST","LGG- SEGA","MB- G4","MB- SHH INF", "MB- SHH CHL AD","MB- G3","EPN- RELA","PIN T-  PB A","SUBEPN- PF",
                        "PTPR- B","SUBEPN- ST","EPN- YAP","EPN- PF A","EPN- PF B","EFT- CIC", "ETMR","CNS NB- FOXR2","LYMPHO","HGNET- BCOR","GBM- MES","GBM- RTK II",
                        "GBM- RTK I","LGG- MYB","CONTR- HEMI","HGNET- MN1","GBM- MID","HMB", "EPN- MPE","SCHW","O IDH","A IDH- HG","MNG","A IDH","LGG- GG","CN",
                        "SUBEPN- SPINE","PIN T-  PB B","PXA","ANA PA","CONTR- INFLAM", "PIN T- PPT","CPH- ADM","CPH- PAP","ENB- A","PITUI","CONTR- PINEAL", "PGG- nC",
                        "LGG- DNT","CHGL","MELAN","PLEX- AD","ENB- B","LIPN", "EPN- SPINE","PTPR- A","SFT HMPC","PLEX- PED A","GBM- RTK III","IHG", "MELCYT","DLGNT",
                        "PITAD- ACTH","PITAD- STH DNS B","PITAD- PRL", "PITAD- FSH LH","PLEX- PED B","EWS","SCHW- MEL","CONTR- ADENOPIT", "LGG- DIG/DIA","PITAD- STH SPA",
                        "PITAD- STH DNS A","CONTR- WM","PLASMA", "CHORDM","RETB","CONTR- PONS","CONTR- HYPTHAL","PITAD- TSH", "bruh"};


#endif