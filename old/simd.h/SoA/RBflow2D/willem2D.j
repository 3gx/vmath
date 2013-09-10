#PBS -S /bin/bash
#PBS -lwalltime=40:05:00
#PBS -lnodes=1:cores8
#PBS -rn
echo Start of Job
export OMP_NUM_THREADS=8
cp -r $HOME/willem2Domp $TMPDIR
(
   cd $TMPDIR/willem2Domp
   ./boutnp
   cp * $HOME/willem2Domp
)
echo End of Job
