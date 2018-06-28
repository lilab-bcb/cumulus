# KCO

[KCO](https://www.broadinstitute.org/klarman-cell-observatory) tools, workflows, and spells

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

#### kco fc-run ####

*kco fc-run* runs a FireCloud method. Features:

1. Uploads local files/directories in your inputs to a Google Cloud bucket updates the file paths to point to the Google Cloud bucket. 
Your sample sheet can point to local file paths and kco run will take care of uploading directories (e.g. fastq directories) and modifying the sample sheet to point to a Google Cloud bucket.
1. Creates or uses an existing workspace.
1. Uses the latest version of a method unless the method version is specified.


##### Options #####
<table>
<tr><td><b>-m<br/>--method</b></td><td>Method namespace/name (e.g. regev/cellranger_mkfastq_count). A version can optionally be specified (e.g. regev/cell_ranger_mkfastq_count/4), otherwise the latest version of the method is used.</td></tr>
<tr><td><b>-w<br/>--workspace</b</td><td>Workspace name (e.g. foo/bar). The workspace will be created if it does not exist</td></tr>
<tr><td><b>-i<br/>--inputs</b></td><td>WDL input JSON. Can be specified as file or string</td></tr>
<tr><td>-n<br/>--no-upload</td><td>Do not upload files/directories to the workspace Google Cloud bucket</td></tr>

</table>

Example: Upload BCL directories and sample sheet, convert sample sheet paths to gs:// URLs, and run cell ranger mkfastq and count.

example_sample_sheet.txt:

<table>
<tr><td>Sample</td><td>Flowcell</td><td>Lane</td><td>Index</td></tr>
<tr><td>s1</td><td>/mylocalpath/flowcell1</td><td>1-3</td><td>index1</td></tr>
<tr><td>s2</td><td>/mylocalpath/flowcell1</td><td>4-6</td><td>index2</td></tr>
<tr><td>s1</td><td>/mylocalpath/flowcell2</td><td>*</td><td>index2</td></tr>
</table>

Note that sample s1 is sequenced on 2 flowcells and all s1 fastq files will be passed to cell ranger count in one command by the pipeline.

inputs.json:

```json
{
"cellranger_count.transcriptome":"GRCh38",
"cellranger_count.sample_sheet":"mylocalpath/example_sample_sheet.txt"
}
```


```
kco run -m regev/cell_ranger_mkfastq_count -i inputs.json -w myworkspace_namespace/myworkspace_name
```


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

