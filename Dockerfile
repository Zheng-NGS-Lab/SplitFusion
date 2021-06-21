FROM python:3.8-slim-buster

LABEL org.label-schema.license="GPL-2.0" \
      org.label-schema.vcs-url="https://github.com/rocker-org/rocker-versioned" \
      org.label-schema.vendor="Rocker Project" \
      maintainer="Carl Boettiger <cboettig@ropensci.org>"

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
	build-essential \
	libc6-dev \
	r-base \
	gcc \
	g++ \
	make \
	bwa \
	bedtools
RUN pip install future \
  && echo "yes" | cpan \
  && Rscript -e "install.packages('data.table')" \
  && Rscript -e "install.packages('plyr')" \
  && perl -MCPAN -e 'install Pod::Usage' \
  && mkdir /annovar && mkdir /SplitFusion && mkdir /data && mkdir /genome
COPY exec /SplitFusion/exec
COPY R /SplitFusion/R
COPY inst /SplitFusion/inst
COPY bin/samtools /usr/bin/samtools
COPY bin/mawk /usr/bin/mawk
WORKDIR /SplitFusion
