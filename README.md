[Read documentation](http://kco-cloud.rtfd.io/)

# KCO

[KCO](https://www.broadinstitute.org/klarman-cell-observatory) tools, workflows, and spells

#### Login to Google ####
Commands that interact with FireCloud require that you've logged in to Google.

If you've done this before you can skip this step - you only need to do this once.

Execute the following command to login to Google.
 ```
gcloud auth login
 ```

Copy and paste the link in your unix terminal into your web browser.
Copy and paste the authorization code in your unix terminal.

#### kco fc_run ####

*kco fc_run* runs a FireCloud method. Features:

1. Uploads local files/directories in your inputs to a Google Cloud bucket updates the file paths to point to the Google Cloud bucket. 
Your sample sheet can point to local file paths and kco run will take care of uploading directories (e.g. fastq directories) and modifying the sample sheet to point to a Google Cloud bucket.
1. Creates or uses an existing workspace.
1. Uses the latest version of a method unless the method version is specified.


##### Options #####
<table>
<tr><td><b>-m<br/>--method</b></td><td>Method namespace/name (e.g. regev/cellranger_mkfastq_count). A version can optionally be specified (e.g. regev/cell_ranger_mkfastq_count/4), otherwise the latest version of the method is used.</td></tr>
<tr><td><b>-w<br/>--workspace</b></td><td>Workspace name (e.g. foo/bar). The workspace will be created if it does not exist</td></tr>
<tr><td><b>-i<br/>--inputs</b></td><td>WDL input JSON. Can be specified as file or string</td></tr>
<tr><td>--bucket-folder</td><td>Store inputs to <folder> under workspace's google bucket</td></tr>
<tr><td>-o<br/>--upload</td><td>Upload files/directories to the workspace Google Cloud bucket and output updated input json (with local path replaced by google bucket urls) to <updated_json>.</td></tr>
</table>

##### Examples #####

Upload BCL directories and sample sheet, convert sample sheet paths to gs:// URLs, and run cell ranger mkfastq and count.

example_sample_sheet.csv:

    ```
    Sample,Reference,Flowcell,Lane,Index,Chemistry
    sample_1,GRCh38,/mylocalpath/flowcell1,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,/mylocalpath/flowcell1,3-4,SI-GA-B8,threeprime
    sample_3,mm10,/mylocalpath/flowcell1,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,/mylocalpath/flowcell1,7-8,SI-GA-D8,fiveprime
    sample_1,GRCh38,/mylocalpath/flowcell2,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,/mylocalpath/flowcell2,3-4,SI-GA-B8,threeprime
    sample_3,mm10,/mylocalpath/flowcell2,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,/mylocalpath/flowcell2,7-8,SI-GA-D8,fiveprime
    ```

Note that sample_1, sample_2, sample_3, and sample_4 are sequenced on 2 flowcells and for each sample, all of its FASTQ files will be passed to cell ranger count in one command by the pipeline.

inputs.json:

```json
{
"cellranger_mkfastq_count.input_csv_file" : "mylocalpath/sample_sheet.csv",
"cellranger_mkfastq_count.cellranger_output_directory" : "gs://url/outputs",
"cellranger_mkfastq_count.delete_input_bcl_directory": true
}
```

```
kco fc_run -m regev/cellranger_mkfastq_count -i inputs.json -w myworkspace_namespace/myworkspace_name --bucket-folder inputs -o inputs_updated.json
```

Upon success, **kco fc_run** returns a url pointing the the submitted FireCloud job. 

If for any reason, your job failed. You could rerun it without uploading files again via the following command:

```
kco fc_run -m regev/cellranger_mkfastq_count -i inputs_updated.json -w myworkspace_namespace/myworkspace_name
```

#### kco sample-sheet ####
*kco sample-sheet* creates a sample sheet.

##### Synopsis #####
```
kco sample-sheet [OPTION]... dir
```

where *dir* can be one or more directories to look for fastq files. Directories are searched recursively for fastq files.

##### Options #####
<table>
<tr><td>-f<br/>--format</td><td>Sample sheet format. Can be directory, r1_r2, r1_r2_i1. For cellranger, use directory, for drop-seq, use r1_r2</td></tr>
<tr><td><b>-o<br/>--output</b></td><td>Output sample sheet</td></tr>
</table>

#### kco fc-upload ####

*kco fc-upload* uploads local files/directories to your FireCloud workspace Google Cloud bucket.
 
Your sample sheet can point to local file paths and kco upload will take care of uploading directories (e.g. fastq directories) and modifying the sample sheet to point to a Google Cloud bucket.


##### Synopsis #####
```
kco fc-upload [OPTION]... file
 ```
where file can be one or more files, such as a sample sheet or an input JSON

##### Options #####
<table>
<tr><td><b>-w<br/>--workspace</b</td><td>Workspace name (e.g. foo/bar). The workspace will be created if it does not exist</td></tr>
<tr><td>--dry-run</td><td>Causes upload to run in "dry run" mode, i.e., just outputting what would be uploaded without actually doing any uploading.</td></tr>

</table>

Example: Upload fastq directories and sample sheet, convert sample sheet paths to gs:// URLs.

example_sample_sheet.txt:

<table>
<tr><td>s1</td><td>/mylocalpath/flowcell1/s1</td></tr>
<tr><td>s2</td><td>/mylocalpath/flowcell1/s2</td></tr>
<tr><td>s1</td><td>/mylocalpath/flowcell2/s1</td></tr>
</table>



```
kco upload -w myworkspace_namespace/myworkspace_name example_sample_sheet.txt 
```

#### kco fc-inputs ####

*kco fc-inputs* generates a stub JSON input file that can be used as input to *kco run*. The JSON file can optionally be based on a published method config.
##### Options #####
<table>
<tr><td><b>-m<br/>--method</b></td><td>Method namespace/name (e.g. regev/cellranger_mkfastq_count). A version can optionally be specified (e.g. regev/cell_ranger_mkfastq_count/4), otherwise the latest version of the method is used.</td></tr>
<tr><td>-c<br/>--config</td><td>Repository config to use for generating input.json stub (e.g. regev/drop-seq-MMUL_8_0_1</td></tr>
<tr><td>-o<br/>--out</td><td>JSON output file</td></tr>
</table>

Example: Generate a stub JSON file based on the published config for running Drop-Seq using the MMUL_8_0_1 genome:

```
kco fc-inputs -m regev/drop-seq -c regev/drop-seq-MMUL_8_0_1
```


