#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt




'/home/sph/code/qe/qe-6.6/bin/pw.x' '-in' 'aiida.in'  > 'aiida.out'
