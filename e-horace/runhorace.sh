#!/bin/bash



mz=90687


for ((i=1; i<=21; ++i))
do


        cat dummy-input | sed -e "s:zmass:$mz:g"\
        		       -e "s:massindex:$i:g"\ > inputfile
        		       
        ./horace < inputfile

       
        mz=$((mz+50))
       
done
exit


