#!/bin/bash

user=mjackson
bindir=/Users/mjackson/Workspace/EIMTomo/Build/Bin
paramdir=/Users/mjackson/Contracts/AFRL-TO81-Extension/TomoData
outputdir=/tmp
cd $bindir



./BasicReconstructionAlgorithm -p ${paramdir}/Params.txt -s ${paramdir}/ScheppLoganSinogram128NoiselessLimitedAngle.bin  -i ${paramdir}/ConstantInitialEst255Slices3Size128X128.bin  -o ${outputdir}/ShepLoganReconstructionSigmaXPoint5P1_2Iter25.bin -l 0.5 -m 1.2 -n 25