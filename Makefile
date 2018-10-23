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

figure-2:
	Rscript figure/2/plot.R > log/figure-2.out 2>&1 &

figure-3:
	Rscript figure/3/plot.R > log/figure-3.out 2>&1 &
