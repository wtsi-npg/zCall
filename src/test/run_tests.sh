#! /usr/bin/env bash

# workaround for annoying NFS error in deleting directories within Python tests

PYTHON_TEST='test/test.py'
DATA_DIR='data'

if [ ! -d $DATA_DIR ]; then
    echo "Output path '$DATA_DIR' does not exist or is not a directory"
    exit 1
elif [ ! -f $PYTHON_TEST ]; then
    echo "$PYTHON_TEST does not exist; should run from src directory of zCall"
    exit 1
else
    python $PYTHON_TEST
    echo "Finished running $PYTHON_TEST";
    rm -Rf $DATA_DIR/output_test_*
fi