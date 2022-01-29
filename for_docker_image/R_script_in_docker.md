#section1. data preparation

what we will need

for the entie procedure we will be needing the following:

1) copy all the .out.dm.ps files in the folder named as [dbCAN_results] generated from the dbCAN annotation platform to the "01_data" folder, in my case the root directory where this "01_data" folder locates is "Desktop/for_docker_image"


#section2. preparation of dockerfile (provided)


#section 3. build and run the image
#modify name of the local root directory if necessary, in my case, the #local root directory is "Desktop/for_docker_image/"
#pwd on mac
#/Users/yuboer/Desktop/for_docker_image

#step1. building image
#note: use the terminal to navigate to the folder (in my case, it is #"Desktop/for_docker_image/" as indicated by the pwd command listed #above) where the Dockerfile is located and build the image with: 

docker build -t cellulosic_genome_annotation .

#step2. running image
#note: 1) using the -v argument tells Docker which local folder to map to #the created folders inside the container
#2) with the command line below, the container will get aceess to data #in the local folder "Desktop/for_docker_image/01_data", and the output #from the workflow "03_output" will be saved locally at 
#"Desktop/for_docker_image/03_output"
docker run -it --rm -v ~/"Desktop/for_docker_image"/01_data:/01_data -v ~/"Desktop/for_docker_image"/03_output:/03_output cellulosic_genome_annotation