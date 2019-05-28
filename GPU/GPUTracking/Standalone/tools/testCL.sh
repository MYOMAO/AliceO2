#!/bin/bash

COMPILER=clang++
LLVM_SPIRV=llvm-spirv

#COMPILER=/usr/lib/llvm/roc-2.1.0/bin/clang++
#COMPILER=/usr/lib/llvm/9/bin/clang++

#COMPILER=/home/qon/llvm/install/bin/clang++
#LLVM_SPIRV=/home/qon/llvm/build-spirv/tools/llvm-spirv/llvm-spirv

INCLUDES="-I ../. -I ../Base -I ../SliceTracker -I ../Common -I ../Merger -I ../TRDTracking -I../ITS -I../dEdx -I$HOME/alice/O2/Detectors/TRD/base/include -I$HOME/alice/O2/Detectors/TRD/base/src -I$HOME/alice/O2/Detectors/ITSMFT/ITS/tracking/include -I$HOME/alice/O2/Common/Constants/include"
DEFINES="-DGPUCA_STANDALONE -DGPUCA_ENABLE_GPU_TRACKER -DGPUCA_GPULIBRARY=OCL -DNDEBUG -D__OPENCLCPP__ -DHAVE_O2HEADERS"

echo Test1 - Preprocess
echo $COMPILER -cl-std=c++ -x cl $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -E ../Base/opencl/GPUReconstructionOCL.cl > test.cl
     $COMPILER -cl-std=c++ -x cl $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -E ../Base/opencl/GPUReconstructionOCL.cl > test.cl
    #Test 1A - Compile Preprocessed
    #$COMPILER -cl-std=c++ -x cl --target=amdgcn-amd-amdhsa -mcpu=gfx906 -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -ferror-limit=1000 -Xclang -finclude-default-header -c test.cl -o test.o
    #exit

echo Test2 - amdgcn
echo $COMPILER -cl-std=c++ -x cl --target=amdgcn-amd-amdhsa -mcpu=gfx906 -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -ferror-limit=1000 -Xclang -finclude-default-header $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -c ../Base/opencl/GPUReconstructionOCL.cl -o test.o
     $COMPILER -cl-std=c++ -x cl --target=amdgcn-amd-amdhsa -mcpu=gfx906 -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -ferror-limit=1000 -Xclang -finclude-default-header $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -c ../Base/opencl/GPUReconstructionOCL.cl -o test.o

echo Test3 - SPIR-V
echo $COMPILER -cl-std=c++ -x cl -emit-llvm --target=spir64-unknown-unknown -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -ferror-limit=1000 -Xclang -finclude-default-header $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -c ../Base/opencl/GPUReconstructionOCL.cl -o test.bc
     $COMPILER -cl-std=c++ -x cl -emit-llvm --target=spir64-unknown-unknown -cl-denorms-are-zero -cl-mad-enable -cl-no-signed-zeros -ferror-limit=1000 -Xclang -finclude-default-header $INCLUDES $DEFINES -Dcl_clang_storage_class_specifiers -c ../Base/opencl/GPUReconstructionOCL.cl -o test.bc
echo $LLVM_SPIRV test.bc
     $LLVM_SPIRV test.bc
