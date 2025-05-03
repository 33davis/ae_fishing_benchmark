#!/bin/bash
#test script

pwd 
#create variable as date set as YYYY_mm_dd
d=$(date +"%Y_%m_%d") 
echo $d

#trying to make a variable directory
mkdir "./$d"

#testing that I can make a sub-directory in variable directory
mkdir "./$d/test_dir"


#testing that I can navigate and appendd
echo "hello world the variable worked!" >  ./$d/test_dir/myfile.txt


#end of script 
