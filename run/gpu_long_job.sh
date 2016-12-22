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
#$ -q gpu_long.q

./../program -r ./../../data/h264_11_00.mkv -in ./../../data/h264_11_18.mkv -type STVSSIM -ffpath ./../../ffmpeg-3.2/ -threads 12
./../program -r ./../../data/h264_11_00.mkv -in ./../../data/h264_11_28.mkv -type STVSSIM -ffpath ./../../ffmpeg-3.2/ -threads 12
./../program -r ./../../data/h264_11_00.mkv -in ./../../data/h264_11_51.mkv -type STVSSIM -ffpath ./../../ffmpeg-3.2/ -threads 12



# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 4
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 5
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 6
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 7
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 8
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 9
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 10
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 11
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 12
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 13
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 14
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 15
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 16
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 17
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 18
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 19
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 20
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 21
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 22
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 23
# ./../program -r ./../../data/h265_11_00.mkv -in ./../../data/h265_11_18.mkv -type PSNR  -ffpath ./../../ffmpeg-3.2/ -threads 24

