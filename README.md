# ClinVar API upload
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) is a public archive of genetic variants and phenotypes, including variant classification with respect to health and disease.

This repository includes scripts to enable easy extraction of classified variants from a [VarSeq Assessment Catalog](https://www.goldenhelix.com/products/VarSeq/) and submission to ClinVar via their [API](https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/).


## Usage

### How to submit variants to ClinVar

#### Workflow overview and background

1. Extract variants from the VarSeq Assessment Catalog

2. Parse, clean and upload variants through the API

> The `txt` extracted from VarSeq is parsed and the required fields extracted to a `json` file needed for the API submission. The `json` will include both novel variants to be submitted and existing variants to be updated.
> 
> Once the variants have been submitted to ClinVar, the corresponding submisison ID and status of the upload are returned by the API, with the list of summary report files saved in a *txt* file. These files are used to annotate the variants in the original *txt* data with SCV accession numbers, in order to retrieve them when variants are to be updated.

3. **IMPORTANT**: Annotate the uploaded variants with SCV accession number and add the annotated variants files to the data folder for the following upload.

> This last step is paramount, as the annotated variants-file serves as a record of variants that have already been submitted to ClinVar. Further, the SCV accession number is required when the variants are to be updated - otherwise the ClinVar API will interpret the updated variant as a new submission.

Disclaimer: all uploaded data are "anonymised" by replacing the "Last Edited" date field in VarSeq with the date of extraction of the catalog.

ClinVar performs a weekly (or otherwise bi-weekly) upload on Sundays, meaning the data uploaded via this tool will only be visible on their website from the following Monday. 

#### 1. Extract variants from the VarSeq Assessment Catalog

##### Obtain VarSeq Assessment Catalog (variant database) URL
* Launch VarSeq (use same version as the interpreters).
* Open a recent project.
* Open Data Source Library. `Menu > Tools > Manage Data Sources > Catalogs`
* Select the catalog and look into the "Information" details view under "Url").


##### Exract variants
Obtain txt file with variants from the VarSeq Catalog
```shell
# Remember to use the most recent varseq version of gautil

# Extract variants from the database using this command
# /path/to/varseq/latest_version/gautil writetsf <txsource_url> output_annotation_file.tsf

# Transform the .tsf to .txt format, i.e.
/path/to/varseq/latest_version/gautil writetxt output_annotation_file.tsf
# Then transfer the resulting .txt file to the local working directory
```

#### 2. Parse and upload variants through the API

##### Clone repository to your local computer

##### Obtain API key
* If your organisation is not registered yet, obtain the key for the API [here](https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/).
* Save the key locally in the cloned repository as "clinvar.key" (see script default = `".\clinvar.key"`).

##### Run `main` script
Use the `main.py` script to parse and upload variants.

The script creates a cleaned version of the catalogue file, `temp\output_annotation_file_cleaned.txt`, removing duplicates and handling variants present on more than one alternate allele (only one allele can be uploaded to ClinVar at a time). The cleaned file is also required later for the annotation of variants in the final step. 

Finally, the list of submission IDs is stored in files as records. When the variants are to be updated, the `check_for_updates.py` script allows to divide the list of variants into those that have already been uploaded, and check if they require to be updated, and the novel variants. These are then uploaded separately.

Output produced by `main.py`
* Successful submission:

    * a `temp` folder containing all the following files;
    * cleaned inputted file of variants;
    * file containing all uploaded, formatted haplotypes;
    * individual *json* files, one per batch, containing the submitted data, saved under the submission ID of each batch, indicating whether the data are variants or haplotypes and if novel or updates (e.g. `SUBid_[variants|haplotypes]_[novel|update].json`);
    * a `summaries_list_[novel|update].txt` file is generated. It contains the paths to the summary report jsons of all batches submitted succesfully. The "novel" file is used to [annotate the data](#annotate-data), following the json retrieval.
* Failed submission: a file containing the errors that led to the failure will be created for each failed batch. In this case, the failed batch is not submitted to the API.

```shell
# See the main script help page
python3 main.py --help

# Example submission with dry run
python3 main.py -f ..\data\output_annotation_file.txt --reference_variants ..\data\variants_uploaded_annotated_{date}.tsv --reference_haplotypes ..\data\haplotypes_uploaded_annotated_{date}.tsv --date_of_extraction {date} -n
```

#### 3. Data annotation

##### Check submission status
The ClinVar API takes a while to process the submitted variants. It is possible to check the status of the individual batch submission using the `check_submission.py` script, providing the associated submission id.

If the submission processing is complete, the script will return a hyperlink to download the summary report json for the specified submission ID. This file **must** be dowloaded and stored in the `temp` directory as it is required for data annotation.

```shell
# See the script help page
python3 check_submission.py --help

# Example use with dry run
python3 check_submission.py -i SUB00000001 -n
```

##### Annotate data

> **N.B.:** The final annotated file must be stored in the data folder in order to be used for the following upload of variants!

Upon successful submission, the variants in the `temp\output_annotation_file_cleaned.txt` variants file must be annotated with the corresponding ClinVar SCV accession numbers. This allows to later update and delete variants via their SCV accession number.

The SCV accession numbers are stored in the summary report json files; in order to annotate the cleaned variants file, the original variants *txt* file is to be passed together with `summaries_list_novel.txt` and the initial reference file from the GitLab repository (stored in `data` folder). This process has to be executed separately for variants and haplotypes with the corresponding references. 

The *--date_of_extraction* parameter must be set to the date of extraction of the variants from the VarSeq database, as this is then stored in the annotated variants file and allows to keep track of the last time the variants were modified in the database.

Steps to follow to ensure correct data annotation following a live run:
* Dowload of summary report json files to `temp` directory (refer to [Check submission status](#check-submission-status));
* Run the `annotate_data.py` script (twice if haplotypes were also uploaded);
* When the files have been annotated, store them (`data\haplotypes_uploaded_annotated_{new-date}.txt`, `data\variants_uploaded_annotated_{new-date}.txt`) in the data folder, replacing the previous reference files.

These steps are essential to keep track of each upload and the variants already present on the online database.

```shell
# See the script help page
python3 annotate_data.py --help

# Example use
# variants
python3 annotate_data.py -s ..\temp\summaries_list_novel.txt -f ..\temp\output_annotation_file_cleaned.txt -r ..\data\variants_uploaded_annotated_{date}.tsv --date_of_extraction {date}
# haplotypes
python3 annotate_data.py -s ..\temp\summaries_list_novel.txt -f ..\temp\haplotypes_uploaded.txt -r ..\data\haplotypes_uploaded_annotated_{date}.tsv -t haplotypes --date_of_extraction {date}
```


### How to delete variants from ClinVar
To delete / remove from ClinVar any variants previously submitted, the `main.py` script can be run with argument `-s delete`. The tab separated file must contain a column, `SCV`, containing the SCV accession number for all variants to be deleted.

```shell
# Example use
python main.py -f ..\data\variants_to_delete.txt -s delete
```


## Questions
All enquiries about the code can be directed to @lprob
