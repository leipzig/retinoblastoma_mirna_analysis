#!/usr/bin/env python
import fileinput

for line in fileinput.input():
    #FGC_36_Gang_Stoff_Laz:1:1:1166:591:GTA:hhh
        parts = line.rstrip().split(':')
        print '@' + ':'.join(parts[0:(len(parts)-3)])
        print parts[len(parts)-2]
        print '+'
        print parts[len(parts)-1]
