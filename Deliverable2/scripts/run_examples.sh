qsub -q short_cpuQ -v MATRIX_FILE="matrices/30k_0p05.mtx" scripts/del2_strong_scaling.pbs

qsub -q short_cpuQ scripts/del2_weak_scaling.pbs