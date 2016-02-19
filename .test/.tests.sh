#! /bin/bash

answerfile='../answers.yml'
test -e $answerfile || { echo "No answers.yml file!";  exit 1; }

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
source ssshtest

set -o nounset

STOP_ON_FAIL=1

answers=(1079 0.384000 850 "PRAME" "chr22:0-16150259")

for answer_num in $(seq 1 5); do
    truth=${answers[$(($answer_num - 1))]}
    run "answer-$answer_num" true
    assert_equal $truth $(grep "answer-$answer_num:" $answerfile \
        | cut -f2 -d' ')
done

