# @ node = 2
# @ tasks_per_node = 32
# @ parallel_threads = 1
# @ task_affinity = core(1) 
# @ notification=never
# @ input=/dev/null
# @ output=RB.$(jobid).out
# @ error=RB.$(jobid).err
# @ wall_clock_limit = 01:00:00
# @ job_type=parallel
# @ network.MPI = sn_all,not_shared,US
# @ queue
export XLSMPOPTS=stack=100000000
cd $HOME/willem
time ./boutnp_24h
