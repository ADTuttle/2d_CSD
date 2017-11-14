#!/bin/bash
#make csd
make debug
./csd -malloc_log -malloc_debug 
#./csd
#./csd -ksp_monitor 
#lldb ./csd
rm csd
