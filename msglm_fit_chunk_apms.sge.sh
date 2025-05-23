
#!/fs/home/stukalov/gentoo/bin/bash
#$ -pe openmp 8
#$ -R y
#$ -l h_rt=2:00:00
#$ -l h_vmem=12G
#$ -N vgirault_vzvapms_msglm
#$ -S /fs/home/stukalov/gentoo/bin/bash
#$ -m ae
#$ -j y
##$ -w v
#$ -o /fs/pool/pool-innate-analysis/scratch/stukalov/logs/vgirault_vzvapms/$JOB_NAME_$JOB_ID_$TASK_ID.log

USER_SCRATCH_PATH=/fs/home/stukalov/scratch
SCRIPTS_PATH=/fs/home/stukalov/projects
PROJECT_ID=vgirault_vzvapms
PROJECT_SCRIPTS_PATH=$SCRIPTS_PATH/adhoc/$PROJECT_ID
# disable MKL threading
export MKL_NUM_THREADS=1
export TMPDIR=$USER_SCRATCH_PATH/Rtmp_qsub
#JOB_ID=1
#JOB_NAME="agebhardt_trachea_msglm"
#SGE_TASK_ID=3
#export NSLOTS=8

EPREFIX="/fs/home/stukalov/gentoo"

# do it!
$EPREFIX/usr/bin/Rscript $PROJECT_SCRIPTS_PATH/msglm_fit_chunk.R \
    $PROJECT_ID $JOB_NAME 20180301 $JOB_ID $SGE_TASK_ID
