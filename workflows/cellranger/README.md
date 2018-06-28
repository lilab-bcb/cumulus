#### CellRanger ####

Broad users can add all tools to their path by executing:
 ```
source FIXME
 ```
 
#### Demultiplex ####
Demultiplex each flowcell using [cellranger mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)



#### Sample Sheet ####
Generate a sample sheet. A sample sheet is a tab separated file with one sample name and fastq directory per line. 
Please note that if same library was sequenced on multiple flowcells, you will have one line for each sample/flowcell.

kco sample_sheet /broad/hptmp/mydir --out sample_sheet.txt


#### Create FireCloud Workspace ####

[Create a FireCloud account](https://software.broadinstitute.org/firecloud/documentation/article?id=6816) if you have not registered for one.
 
[Create a FireCloud workspace]((https://software.broadinstitute.org/firecloud/documentation/article?id=10746)) or use an existing workspace.



#### Login to Google ####
Authenticate with Google. If you've done this before you can skip this step - you only need to do this once.


Execute the following command to login to Google.
 ```
gcloud auth login
 ```

Copy and paste the link in your unix terminal into your web browser.
Copy and paste the authorization code in your unix terminal.


#### Upload ####
Upload the sample sheet and fastq directories to your FireCloud workspace bucket.

kco upload sample_sheet.txt -w workspace_namespace/workspace_name

#### Configure Method ####





#### Launch Count ####
