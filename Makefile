SLURM40=
# uncomment to submit to slurm
# adjust the number of cores in the scripts accordingly
# SLURM40=srun -c 40 --mem=40GB


simulate-time-point-sets:
	${SLURM4} Rscript simulations/contribution-of-points.R	> log/$@.out 2>&1 &

fit-experiment:
	${SLURM40} Rscript  experiment/fit.R  > log/$@.out 2>&1 &

compute-conf-int:
	${SLURM40} Rscript  experiment/ci.R > log/ci.out 2>&1 &

prepare-slamseq-data:
	 Rscript  slamseq/01-load.r > log/ci.out 2>&1 &

fit-slamseq:
	 Rscript slamseq/02-fit.R  > log/$@.out 2>&1
	 Rscript slamseq/03-make-csv.r  > log/$@.out 2>&1

figure-2:
	Rscript figure/2/plot.R > log/figure-2.out 2>&1 &

figure-3-4:
	Rscript figure/slamseq/plot.R > log/figure-3-4.out 2>&1 &

figure-6:
	Rscript figure/conventional/plot.R > log/figure-6.out 2>&1 &
