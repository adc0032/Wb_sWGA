#!/usr/bin/env python
import os
import sys
import subprocess

from bx.align import maf


def main(in_file, org, chrom, start, end):
    index_file = in_file + ".index"
    if not os.path.exists(index_file):
        build_index(in_file)

    region_name = "%s.%s" % (org, chrom)
    start = int(start)
    end = int(end)
    index = maf.Indexed(in_file, index_file)
    for align in index.get(region_name, start, end):
        region_align = align.slice_by_component(region_name, start, end)
        seqs_by_org = dict()
        for component in region_align.components:
            seqs_by_org[component.src] = component.text
        print seqs_by_org

main("bmal_wb.maf", "Wb_PNG_Genome_assembly_pt22", "PairedContig_6005", 86, 87)


def build_index(in_file):
    cl = ["maf_build_index.py", in_file]
    subprocess.check_call(cl)

if __name__ == "__main__":
    main(*sys.argv[1:])
