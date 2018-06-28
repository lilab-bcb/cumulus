#### CellRanger ####

Broad users can add necessary command line tools to their path by executing:
 ```
source /seq/regev_genome_portal/conda_env/kco_env/kco.sh
 ```
 
#### Demultiplex ####
Demultiplex each flowcell using [cellranger mkfastq](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq)


#### Create Sample Sheet ####
Generate a sample sheet. A sample sheet is a tab separated file with one sample name and fastq directory per line. 
Please note that if same library was sequenced on multiple flowcells, you will have one line for each sample/flowcell.

kco sample_sheet /broad/hptmp/mydir1 /broad/hptmp/mydir2  --out sample_sheet.txt

The following example sample sheet indicates that sample s1 was sequenced on two flowcells, while s2 was sequenced on one flowcell:
<table>
<tr><td>s1</td>/broad/hptmp/mydir1/s1</tr>
<tr><td>s2</td>/broad/hptmp/mydir1/s2</tr>
<tr><td>s1</td>/broad/hptmp/mydir2/s1</tr>
</table>

#### Create FireCloud Workspace ####

[Create a FireCloud account](https://software.broadinstitute.org/firecloud/documentation/article?id=6816) if you have not registered for one.
 
[Create a FireCloud workspace]((https://software.broadinstitute.org/firecloud/documentation/article?id=10746)) or use an existing workspace.

Make note of your workspace id. For example, in the following screenshot, help-gatk/Pre-processing_hg38_v2 is the workspace id.

![Workspace Screenshot](https://klarman-cell-observatory.github.io/KCO/workflows/cellranger/images/workspace.png)

#### Login to Google ####
If you've done this before you can skip this step - you only need to do this once.

Execute the following command to login to Google.
 ```
gcloud auth login
 ```

Copy and paste the link in your unix terminal into your web browser.
Copy and paste the authorization code in your unix terminal.


#### Upload Files ####
Upload the sample sheet and fastq directories to your FireCloud workspace bucket.

kco fc_upload sample_sheet.txt -w help-gatk/Pre-processing_hg38_v2

*fc_upload* will automatically upload all fastq directories referenced in your sample sheet and will convert the paths to the fastq directories to Google Cloud bucket URLs.


#### Configure FireCloud Method ####

Click on "Method Configuration" tab. Then click on "Import Configuration.
Click on "Import from Method Repository" in the dialog box.
Under Public Methods, use the search box to find regev/cellranger_count
Click "Use Blank Configuration" and then click "Import Method"
Click "Edit Configuration" and scroll down to the Inputs section
Enter the gs:// URL to your sample sheet. Please note that you must put quotes around the URL, e.g. "gs://fc-fb46edeb-f05e-45c9-88f7-9ef83d8f1aac/sample_sheet.txt"
Change any additional parameters.
Click "Save".

#### Launch FireCloud Method ####
Under the "Method Configurations" tab you will find the Launch Analysis button.
