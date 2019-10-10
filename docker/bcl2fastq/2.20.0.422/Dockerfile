FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

# See license at https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf

RUN apt-get update && \
apt-get install --no-install-recommends -y build-essential dpkg-dev gnupg lsb-release procps && \
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
pip install pandas

RUN mkdir /software
ADD bcl2fastq2_0v2.20.0.422-2_amd64.deb /software

RUN dpkg -i /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb && rm /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb
ENV PATH=/software:$PATH
