#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt


'/home/sphuber/code/qe/qe-6.1/bin/pw.x' '-in' 'aiida.in'  > 'aiida.out' 
