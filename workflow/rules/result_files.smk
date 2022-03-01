# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Jonas A"
__copyright__ = "Copyright 2021, Jonas A"
__email__ = "jonas.almlof@igp.uu.se"
__license__ = "GPL-3"


rule copy_results_files:
    input:
        input_files,
    output:
        output_files,
    run:
        import subprocess
        i = 0
        for file in input[0]:
            subprocess.run(["cp", file, output[0][i]])
            i += 1
