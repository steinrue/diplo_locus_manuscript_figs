#!/bin/bash

#Request 1 processor
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

# Request 1 hour of wall clock time
#SBATCH --time=48:00:00

# Request memory per process
#SBATCH --mem=8gb

# Job Array###SBATCH --array=3-5

# Request regular output and terminal output go to the same file
#SBATCH --output=%x.o%j-%a
#SBATCH --error=%x.o%j-%a


# some sanity check
echo $( date ) 'Starting...'


#### Pipeline to perform simulation for all replicate in a given scenario, ####
### run DL on each rep, meanwhile collect data from all selected sites for ###
### each rep, and run DiploLocus on it.                                        ####


# path='/gpfs/data/steinruecken-lab/XiaohengStuff/Diffusion_spect/Simulations/singleLocus/'
path=$( pwd )

slimWrapper=${path}"run_MSinit_slim_reps_SV.py"

parsingScript=${path}'parseSlim4_TrackSel_var-t.py'

cd $path

## variables commented out below will be passed from sbatch commandline
lambda=1
sampledTimes=9
sampSize=40

## for testing
## below are variables needed from the command line
# frac_s=$1
# h=$2
# part=40
# start=$3
# end=$4

s='0'${frac_s}
if [[ $h == 0.5 ]]; then
	cond='s'${frac_s}'_h.5'
else
	cond='s'${frac_s}'_h'${h}
fi

# some sanity checks
echo 'frac_s='${frac_s}', s='${s}', h='${h}', cond='${cond}
echo "job array id="${SLURM_ARRAY_JOB_ID}
echo "trying out other variable name: \$SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID}
start=$( expr ${SLURM_ARRAY_TASK_ID} \* $part )
end=$( expr ${SLURM_ARRAY_TASK_ID} \* $part + $part - 1 ) 

# get seeds
seedfile=${path}'slim4_input/200kb_HC_'${cond}'_ConstN_t'${sampledTimes}'x500gen_simSeeds.txt'
# # get the planned total # of reps:
# Nrep=$( cat ${seedfile} | awk '($1 ~ /[0-9]+/)' | wc -l  )
# just set it to 200 for now
Nrep=200

date
echo 'From rep '${start}' to rep '${end}'. '${Nrep}' reps in total'
echo ""

# neut & selection use the same script
slimin=${path}'Slim4_200kb_HC_sel-noTrack_stdVarSFS_ConstN1e4_t9x500gen_var-n_Lambda1.eidos'
# no need to extract seed
## just to make sure other array jobs "wait" for the first one a little bit
## for the directories to be made; seed files already exist
if [[ $start -gt 0 ]]; then
  sleep 5
fi

# run msprime and slim in a batch (save seed-reading time)
python $slimWrapper $slimin $cond $start $end "MAF.05" "reparse"

for rep in $( seq $start $end ); do
  countfile=${path}'200kb_HC_'${cond}'_ConstN_t'${sampledTimes}'x500gen/count/HC_'${cond}'_ConstN_t'${sampledTimes}'x500gen_n'${sampSize}'_rep'${rep}'.count'
	if [[ ! -e $countfile ]]; then
	  python $slimWrapper $slimin $cond $rep $rep "MAF.05" "reparse"
	fi
	# now let's run stuff
	outprefix=${path}'200kb_HC_'${cond}'_ConstN_t'${sampledTimes}'x500gen/likelihood/HC_'${cond}'_ConstN_t'${sampledTimes}'x500gen_n'${sampSize}'_rep'${rep}'_DL'

	# do geom grid
	gridInfo='51xgeom75e-2'
	# pwd
	if [[ ! -f ${outprefix}'_MAF.05_'${gridInfo}'_Unif_off-grid_maxLLs.txt' ]]; then
		# record the command
		date
		echo DiploLocus likelihood -i $countfile --u01 1.25e-8 --Ne 1e4 \
				--sample_times 0,500,1000,1500,2000,2500,3000,3500,4000 --minMAF 0.05 \
				--geom_s2_range="-0.75,0.75,50" --fix_h $h --init uniform --gzip_surface \
				--get_off_grid_max --get_MLR --get_chi2_pval -o ${outprefix}'_MAF.05_'${gridInfo}'_Unif'
		echo ""
		DiploLocus likelihood -i $countfile --u01 1.25e-8 --Ne 1e4 \
				--sample_times 0,500,1000,1500,2000,2500,3000,3500,4000 --minMAF 0.05  \
				--geom_s2_range="-0.75,0.75,50" --fix_h $h --init uniform --gzip_surface \\\
				--get_off_grid_max --get_MLR --get_chi2_pval -o ${outprefix}'_MAF.05_'${gridInfo}'_Unif'
	fi
	echo ""
# end of looping through reps
done

date
echo "Done."

sacct -j $SLURM_JOB_ID -o JobID,JobIDRaw,JobName,AveCPU,ntasks,AllocCPUS,Elapsed,State,ExitCode,ReqMem,ConsumedEnergy,MaxDiskRead,MaxDiskWrite
#Partition,MinCPU,MinCPUTask,MaxVMSize,AveVMSize,MaxRSS,AveRSS,,AveCPUFreq,ReqCPUFreqMin,ReqCPUFreqMax,AveDiskRead,AveDiskWrite