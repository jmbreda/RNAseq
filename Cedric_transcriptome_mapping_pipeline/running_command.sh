#!/bin/sh
snakemake  -s Snakefile  -j 80 --restart-times 2 --cluster-config cluster.json --ri --cluster "sbatch --cpus-per-task {cluster.n}  --time {cluster.time} --mem {cluster.mem} --qos serial"
#snakemake  -s Snakefile  -j 999 --restart-times 3 --cluster-config cluster.json --use-envmodules --cluster "sbatch --cpus-per-task {cluster.n}  --time {cluster.time} --mem {cluster.mem} --qos serial"


