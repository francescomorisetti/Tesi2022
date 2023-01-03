#!/bin/bash


decdummy=107
mz=91.$decdummy


for ((i=1; i<=15; ++i))
do


        cat dummy-input | sed -e "s:zmass:$mz:g"\
        		       -e "s:massindex:$i:g"\ > inputfile
        		       
        ./horace < inputfile

        decdummy=$((decdummy+10))
        mz=91.$decdummy
       
done
exit


