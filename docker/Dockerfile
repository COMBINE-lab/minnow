# image: COMBINE-lab/minnow
# This dockerfile is based on the one created by
# Titus Brown (available at https://github.com/ctb/2015-docker-building/blob/master/salmon/Dockerfile)
#FROM insilicodb/bioconductor 
FROM ubuntu:16.04

ENV PACKAGES git gcc make g++ cmake  liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev curl unzip autoconf vim wget time bzip2
    ENV MINNOW_VERSION 1.0

# TwoPaCo
RUN wget https://github.com/01org/tbb/releases/download/2018_U4/tbb2018_20180411oss_lin.tgz && \
tar zxf tbb2018_20180411oss_lin.tgz && \
rm tbb2018_20180411oss_lin.tgz && \
cd tbb2018_20180411oss && \
cp -r include /usr/local && \
cp -r lib/intel64/gcc4.7/* /usr/local/lib && \
cd .. && \
rm -rf tbb2018_20180411oss

RUN git clone https://github.com/fataltes/TwoPaCo.git

RUN cd TwoPaCo && \
	mkdir build && \
	cd build && \
	cmake ../src && \
	make

# Minnow
RUN git clone https://github.com/COMBINE-lab/minnow.git
RUN cd minnow && \	
	mkdir build && \
	cd build && \
	cmake ../ && \
	make

RUN ldconfig
