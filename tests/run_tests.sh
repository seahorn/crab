#!/bin/bash  

usage() {
    echo "Usage: $0 EXPECTED_OUTPUT BUILD_DIR"
    echo "where"
    echo "   EXPECTED_OUTPUT = expected_results.[apron|pplite|elina|boxes].out"
    echo "   BUILD_DIR is the build directory (e.g., build) where tests are located"
}

if [ $# -ne 2 ]; then
    echo "ERROR: $0 expects 2 parameters but given $#"
    usage
    exit 1
fi

OLDLOG=$1
DIR=$2

if [ ! -d "$DIR/test-bin" ]; then
    echo "$DIR/test-bin does not exist"
    exit 1
fi

if [ ! -f "$OLDLOG" ]; then
    echo "$OLDLOG does not exist"
    exit 1
fi

## From OLDLOG we can extract which external library tests (if any)
## should be included
lib=$(basename -- "$OLDLOG")
lib="${lib##expected_results.}"
lib="${lib%.out}"

if [[ $lib != "apron" ]] && [[ $lib != "pplite" ]] && [[ $lib != "boxes" ]] && [[ $lib != "elina" ]] && [[ $lib != "out" ]]; then
    echo "First parameter should be expected_results.[apron|pplite|elina|boxes].out but it is $OLDLOG"
    exit 1
fi    


#DIFF=/Applications/DiffMerge.app/Contents/MacOS/DiffMerge
DIFF=diff

timestamp=$(date +"%m_%d_%y.%H_%M")  
NEWLOG=results_${timestamp}.out


## Run all the tests
for test in $DIR/test-bin/*
do
  if [[ $lib == "out" ]]; then
      if [[ $test != *"apron"* ]] && [[ $test != *"pplite"* ]] && [[ $test != *"boxes"* ]] && [[ $test != *"elina"* ]]; then    
	  echo "Running $test"
	  echo "=== Begin $test ===" >> $NEWLOG
	  $test >> $NEWLOG 2>/dev/null
	  echo "=== End $test ===" >> $NEWLOG
      fi
  else
      # lib should be apron, pplite, boxes, or elina
      if [[ $test == *"$lib"* ]]; then    
	  echo "Running $test"
	  echo "=== Begin $test ===" >> $NEWLOG
	  $test >> $NEWLOG 2>/dev/null
	  echo "=== End $test ===" >> $NEWLOG
      fi
  fi
done

LOGDIR=$(mktemp -d "${TMPDIR:-/tmp/}$(basename $0).XXXXXXXXXXXX")
LOG="$LOGDIR/log.txt"

## Diff the output of the tests with the expected output
$DIFF --suppress-common-lines \
      --ignore-matching-lines="CRAB WARNING:*" \
      --ignore-matching-lines="=== *" \
      $OLDLOG $NEWLOG >& $LOG
STATUS=$?

###################################################################
## Comment this line if we want to keep the generated output
###################################################################
rm -f $NEWLOG

if [ $STATUS -eq 0 ]; then
    echo "All tests passed successfully!"
    exit 0
else
    echo "Some tests produced unexpected output:"
    cat $LOG
    exit 1  
fi
