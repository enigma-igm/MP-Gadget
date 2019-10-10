#! /bin/bash

#
# this script runs a benchmark for the code rev directory
# generated by build-MPGadget.sh
# must be run from the dir where the code is. link it there!


code=$1
suite=$2
jobscript=$3

if [ "x$2" == "x" ]; then
    echo usage $0 code suite
    echo example: $0 coriknl-xxxxxxxx dm-50-512 job.coriknl
    echo suite must be declared in benchmarks directory and code must have
    echo been created by build-MPGadget.sh
fi

# make sure there is no tailing /
code=${code///}

trap exit ERR

prefix=$CSCRATCH/mp-gadget-benchmarks/$code
codedir=`readlink -f $code`

builddir=$codedir/build

logdir=$PWD/logs

mkdir -p $logdir

(cd $codedir/benchmarks; bash gen-bench.sh $suite $prefix/$suite)

export prefix suite logdir code jobscript

cat $prefix/$suite/$jobscript | sbatch -A m3058 --qos premium $C \
    -o $logdir/${code}-${suite}-${jobscript}.o%j