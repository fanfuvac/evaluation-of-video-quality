#!/bin/sh

#  ===========================================================================
# |                                                                           |
# |             COMMAND FILE FOR SUBMITTING SGE JOBS                          |
# |                                                                           |
# |                                                                           |
# | SGE keyword statements begin with #$                                      |
# |                                                                           |
# | Comments begin with #                                                     |
# | Any line whose first non-blank character is a pound sign (#)              |
# | and is not a SGE keyword statement is regarded as a comment.              |
#  ===========================================================================

# Request Bourne shell as shell for job
#$ -S /bin/sh

# Execute the job from the current working directory.
#$ -cwd

# Defines  or  redefines  the  path used for the standard error stream of the job.
#$ -e .

# The path used for the standard output stream of the job.
#$ -o .

# Do not change.
#$ -pe ompi 1

# Do not change.
#$ -q gpu_01.q


./../program -r ./../../data/h265_fHD_crf22.mkv -in ./../../data/h264_fHD_crf26.mkv -type STVSSIM CUDA -ffpath ./../../ffmpeg-3.2/ -CUDA_threads 64

