#!/usr/bin/python

import sys
import os

if len(sys.argv) != 2:
    print sys.argv[0], "file.txt"
    exit(1)

fin = open(sys.argv[1], "r")

basename, ext = os.path.splitext(os.path.basename(sys.argv[1]))

fr3_out  = open(basename + "_fr3.bed", "w")
cdr3_out = open(basename + "_cdr3.bed", "w")

sequence = ""
start = False
for line in fin.readlines():
    if line[0:6] == "Query=":
        start = True

        fr3_start  = -1
        fr3_end    = -1
        cdr3_start = -1
        cdr3_end   = -1

        sequence = line[7:].strip()
    elif not start:
        continue
    else:
        # we are inside a section, let's find the CDR3 coordinates
        fields = line.split()

        if len(fields) < 5:
            continue

        if fields[0] == "CDR3-IMGT":
            cdr3 = True

            cdr3_start = int(fields[2]) - 1
            cdr3_end   = int(fields[3])

            # we only save the FR3 sequences which are followed by a CDR3
            if fr3_start > -1 and fr3_end > -1:
                print >>fr3_out, "%s\t%s\t%s" % (sequence, fr3_start, fr3_end)

            print >>cdr3_out, "%s\t%s\t%s" % (sequence, cdr3_start, cdr3_end)

            start = False # skip next lines until next section

            fr3_start  = -1
            fr3_end    = -1
            cdr3_start = -1
            cdr3_end   = -1

        elif fields[0] == "FR3-IMGT":
            fr3_start = int(fields[1]) - 1
            fr3_end   = int(fields[2])

fin.close()
fr3_out.close()
cdr3_out.close()
