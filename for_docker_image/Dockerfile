#Base image
FROM rocker/tidyverse
#library installation
RUN install2.r --error --deps TRUE dplyr genoPlotR

#create directories
RUN mkdir -p /01_data
RUN mkdir -p /02_code
RUN mkdir -p /03_output

#copy files
#COPY /01_data/ /01_data/
COPY /02_code/docker_MIMAG_annotation_final.R /02_code/docker_MIMAG_annotation_final.R


#run the script
CMD Rscript /02_code/docker_MIMAG_annotation_final.R







