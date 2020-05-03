#!/usr/bin/bash

samplepath="/projectnb/bf528/users/group4/project4/code/salmon-latest_linux_x86_64/bin"
.$samplepath/salmon index -i index -k 31 --gencode -p 4 -t gencode.v33.transcripts.fa.gz
