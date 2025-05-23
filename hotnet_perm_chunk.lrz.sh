#!/bin/bash
#SBATCH -o /gpfs/scratch/pn69ha/%u/logs/vgirault_vzvapms/%x_%j.log
#SBATCH -J vzvapms_hotnet_perm
#SBATCH --get-user-env

#SBATCH --clusters=cm2_tiny
#SBATCH --nodes=4

##SBATCH --clusters=cm2
##SBATCH --partition=cm2_std
##SBATCH --qos=cm2_std
##SBATCH --nodes=3-24

##SBATCH --clusters=cm2
##SBATCH --partition=cm2_large
##SBATCH --qos=cm2_large
##SBATCH --nodes=25-32

#SBATCH --mincpus=2
##SBATCH --mem-per-cpu=2GB
##SBATCH --exclusive=user
#SBATCH --cpus-per-task=10
##SBATCH --ntasks-per-socket=1
#SBATCH --mail-type=end
#SBATCH --mail-user=alexey.stukalov@tum.de
#SBATCH --export=NONE
#SBATCH --time=48:00:00

IMAGES_PATH=$SCRATCH/docker4muc
CHUDIS_PATH=$HOME/projects/cool_chunk_dispatcher

PROJECT_ID=vgirault_vzvapms
HOTNET_VERSION=20210707
CHUDIS_JOBID=${PROJECT_ID}_${HOTNET_VERSION}_perm
NTASKS_PER_NODE=2

module load slurm_setup
module load charliecloud

NQUEUES=$(($SLURM_JOB_NUM_NODES * $NTASKS_PER_NODE))
echo "$NQUEUES chunk dispatcher queue(s) will be started"

for i in $(seq $NQUEUES); do
echo "Starting queue #${i}..."
srun -c10 -n1 --nodes=1 --ntasks-per-node=$NTASKS_PER_NODE --mem-per-cpu=2GB --wait=0 \
     --no-kill --distribution=block --exclusive=user \
     -o $SCRATCH/logs/$PROJECT_ID/${CHUDIS_JOBID}_${SLURM_JOB_ID}_$i.log \
${CHUDIS_PATH}/process_chunks.sh $CHUDIS_JOBID $USER ${SLURM_JOB_ID}_$i \
"ch-run $IMAGES_PATH/archpc.julia \
     -t --no-home --unset-env='*PATH' \
     --set-env=$HOME/projects/adhoc/$PROJECT_ID/hotnet_julia.lrz.env \
     -b $HOME/projects/adhoc:/projects/adhoc \
     -b $HOME/scratch:/scratch -- \
     julia --project=/projects/adhoc/$PROJECT_ID \
     /projects/adhoc/$PROJECT_ID/hotnet_perm_chunk.jl \
     $PROJECT_ID $SLURM_JOB_NAME $HOTNET_VERSION $CHUDIS_JOBID 10" &
done

wait
echo "All queues of job #$CHUDIS_JOBID finished"
