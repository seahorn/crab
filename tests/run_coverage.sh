#!/bin/bash  

if [ $# -ne 2 ]; then
    echo "Usage: $0 BUILD_DIR CRAB_ROOT"
    echo "BUILD_DIR is the build directory (e.g., build) where tests are located"
    exit 1
fi

BUILD_DIR=$1
CRAB_ROOT=$2

if [ ! -d "$BUILD_DIR/test-bin" ]; then
    echo "$BUILD_DIR/test-bin does not exist"
    exit 1
fi

OUTFILE=all.info
REPORT=lcov_report

LCOV=lcov
GENHTML=genhtml

## Run all the tests
echo 'Running tests ...'
for test in $BUILD_DIR/test-bin/*
do
  $test >& /dev/null
done
echo 'done!'

# Generate coverage information
if [ ! -d $BUILD_DIR/lcov_files ]; then
    mkdir $BUILD_DIR/lcov_files
fi    
echo 'Generating coverage information ...'
PATTERNS=
for test in $BUILD_DIR/test-bin/*
do
    test=${test##*/}
    $LCOV -c -d $BUILD_DIR/tests/CMakeFiles/$test.dir/ -b $BUILD_DIR \
	  -o $BUILD_DIR/lcov_files/coverage.$test.info
    $LCOV -e $BUILD_DIR/lcov_files/coverage.$test.info \
	  "${CRAB_ROOT}/include/crab/*" \
	  "${CRAB_ROOT}/lib/*" \
	  -o $BUILD_DIR/lcov_files/$test.info
    cat $BUILD_DIR/lcov_files/coverage.$test.info >> $BUILD_DIR/${OUTFILE}
done
echo 'done!'

echo 'Generating report ...'
$GENHTML  $BUILD_DIR/${OUTFILE} --output-directory $BUILD_DIR/$REPORT
