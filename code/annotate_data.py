import datetime
import json
import csv
import argparse
from pathlib import Path
import datetime
import re
import sys
import helper_functions

def parse_arguments():
    parser = argparse.ArgumentParser(description="Provide cleaned data file and summary json for variant annotation.")
    parser.add_argument(
        "-s",
        "--summary_files",
        help='Path to file containing list of summary report json paths (created from the "main.py" script)',
        type=Path,
        required=True
    )
    parser.add_argument(
        "-f",
        "--input_file",
        help="Path to cleaned file to be annotated with SCV numbers",
        type=Path,
        required=True
    )
    parser.add_argument(
        "-o",
        "--output_file",
        help="Path for outputted file, annotated with SCV numbers",
        type=Path,
        default= None
    )
    parser.add_argument(
        "-r",
        "--reference_file",
        help="Path for reference file (GitLab), annotated with SCV numbers from previous uploads",
        type=Path
        # required=True
    )
    parser.add_argument(
        "-t",
        "--variant_type",
        help="Type of variant to be annotated (default = \"variants\")",
        choices=["variants", "haplotypes"],
        default="variants"
    )
    parser.add_argument(
        "--date_of_extraction",
        help="Date of extraction of variants from VarSeq (API upload date does not reflect this as always done with some delay)",
        required = True
    )
    parser.add_argument(
        "--somatic",
        help="If somatic catalogue is to be annotated, pass this flag",
        action="store_true"
    )
    return parser


def get_new_annotation(summary_report, annotation_dict, date_string):
    """Function to get the SCV accession number for all uploaded variants from the summary report json.
    If the variant was not uploaded but is already present and was not previously annotated, it will
    be this time round and will need to be updated.

    Parameters:
        summary_report: path to the summary report json file returned from "check_submission.py" script.
        annotation_dict: initial dictionary containing annotation from reference file.
        date_string: date of upload, for further annotation.
    
    Return:
        annotation_dict: updated dictionary containing for each submitted variant the corresponding SCV accession number from the current report summary json file parsed.
    """
    with open(summary_report) as submission_response:  # json to be read
        data = json.load(submission_response)
        print(f'Annotating from date: {data["submissionDate"]}, submission: {summary_report}')
    for submission in data["submissions"]:
        try:
            ID, SCV = (
                submission["identifiers"]["clinvarLocalKey"].split("|")[0],
                submission["identifiers"]["clinvarAccession"],
            )  # successful uploads
        except:
            ID = submission["identifiers"]["clinvarLocalKey"].split("|")[0]
            if submission["errors"][0]["output"]["errors"][0][
                "userMessage"
            ].startswith(
                "This record is submitted as novel but"
            ):  # variants present, to be updated
                message = submission["errors"][0]["output"]["errors"][0][
                    "userMessage"
                ].split()
                SCV = message[23]
            else: #other reason for unseccesful upload; variant to be re-submitted
                continue
        annotation_dict[ID] = {"SCV": SCV, "Last Edited": date_string}
    return annotation_dict


def annotate_from_ref(reference_file):
    """Function to extract the SCV accession numbers from the reference file.

    Parameters:
        reference_file: reference file with previously uploaded variants, containing SCV field.
    
    Return:
        partial_annotation: dictionary where variant's ID contains its SCV accession number.
    """
    partial_annotation = {}
    if reference_file:
        with open(reference_file, "r", encoding="utf-8") as ref:
            for line in csv.DictReader(ref, delimiter="\t"):
                if len(line["hgvs c."]) > 0:
                    ID = line["hgvs c."].split(",")[0] #if multiple hgvs c. are present, the first is considered the correct one by default and matching to the reference file
                else:
                    ID = get_id(line)
                if re.search(r"\(\d+\)", ID):
                    ID = helper_functions.fix_hgvs(ID, line["Ref/Alt"].split("/")[1])
                SCV = line["SCV"]
                partial_annotation[ID] = {"SCV": SCV, "Last Edited": line["Last Edited"]}
                # partial_annotation[ID] = {"SCV": SCV, "Last Edited": date_string}
    return partial_annotation


def get_id(variant):
    """If HGVS c. is not present for the variant, the function extracts the variant information
    in the format required by the API to upload it and, consequently, to compare to the
    variants in the reference file.
    
    Parameters:
        variant: variant missing HGVS c. with information to be extracted.
    
    Return:
        variant_check: variant id in correct format.
        """
    ref, alt = variant["Ref/Alt"].split("/")
    variant_check = "_".join(
        [
            f'Chr.{variant["#Chromosome"]}',
            str(
                int(variant["Start"]) + 1
            ),  # in summary report start position is increased by 1
            variant["Stop"],
            ref,
            alt,
        ]
    ) #constructed ID to access dictionary with SCV numbers
    if "-" in variant["Ref/Alt"].split("/")[0]:
        variant_check = "_".join([variant_check, "Insertion"])
    elif "-" in variant["Ref/Alt"].split("/")[1]:
        variant_check = "_".join([variant_check, "Deletion"])
    return variant_check


def get_summary_filepaths(summaries_list, variant_type = "variants"):
    """Function to extract the file names for all submitted batches' summaries.
    
    Parameter:
        summaries_list: txt file containing the list of filepaths of submissions' summary jsons.
        variant_type: identification of type of variant to be annotated, either "variants" of "haplotypes"; used to parse the summaries_list file.
    
    Return:
        filenames: list of file paths to summary reports.
    """
    with open(summaries_list) as f:
        lines = f.readlines()
    filenames = []
    for line in lines:
        if variant_type in line:
            summary_report = line.split("\n")[0].split(" ")[0]
            filenames.append(summary_report)
    return filenames


def annotate_file(summary_reports, input_file, reference_file, output_file, variant_type, date_string):
    """Function to annotate the inputted, cleaned file with SCV accession numbers from upload to ClinVar.
    The function writes a new file in the same tab-delimited format with the addition of "SCV" column.

    Parameters:
        summary_report: summary report json file to access SCV numbers.
        input_file: original file of input (cleaned).
        output_file: final filename with annotated variants.
        variant_type: whether to annotate 'variants' or 'haplotypes'.
        date_string: date of extraction of variants from VarSeq.
    """
    annotation = annotate_from_ref(reference_file)
    for f in get_summary_filepaths(summary_reports, variant_type):
        annotation = get_new_annotation(f, annotation, date_string)
    all_rows = []
    if variant_type == "haplotypes": #read from the old annotated file the previous uploaded haplotypes to then be copied in the new file
        if reference_file:
            with open(reference_file, encoding="utf-8", mode="r+") as old_haplotypes:
                csv_reader = csv.DictReader(old_haplotypes, delimiter="\t")
                for row in csv_reader:
                    all_rows.append(row)
    with open(input_file, encoding="utf-8", mode="r+") as non_annotated:
        csv_reader = csv.DictReader(non_annotated, delimiter="\t")
        if "SCV" not in csv_reader.fieldnames:
            csv_reader.fieldnames.append("SCV")
            field_names = csv_reader.fieldnames
        for row in csv_reader:
            if row["hgvs c."]: #if hgvs c. field present for the variant this is used to find corresponding SCV number
                hgvs_field = row["hgvs c."].split(",")[0]
                hgvs_field_new = hgvs_field
                if re.search(r"\(\d+\)", hgvs_field):
                    hgvs_field_new = helper_functions.fix_hgvs(hgvs_field, row["Ref/Alt"].split("/")[1])
                    updated_hgvs = True
                if (hgvs_field in annotation.keys()) or (hgvs_field_new in annotation.keys()):
                    row["SCV"] = annotation[hgvs_field]["SCV"]
                    row["Last Edited"] = annotation[hgvs_field]["Last Edited"]
                    row["hgvs c."] = hgvs_field
            else: #otherwise the ID has to be constructed with the variant's information
                variant_check = get_id(row)
                if variant_check in annotation.keys(): #check whether variant was uploaded/updated
                    row["SCV"] = annotation[variant_check]["SCV"]
                    row["Last Edited"] = annotation[variant_check]["Last Edited"]
            all_rows.append(row)
            
    with open(output_file, encoding="utf-8", mode="w", newline="") as outfile:
        csv_writer = csv.DictWriter(
            outfile, fieldnames=field_names, delimiter="\t"
        )
        csv_writer.writeheader()
        csv_writer.writerows(all_rows)


if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()
    if len(args.date_of_extraction) == 10:
        try:
            datetime.datetime.strptime(args.date_of_extraction, '%Y-%m-%d')
            # date_check = bool(date_parser.parse(args.date_of_extraction))
            date_of_interest = args.date_of_extraction
        except ValueError:
            print("Error: Date of extraction must be in format %YYYY-%mm-%dd")
            sys.exit(1)
    else:
        print("Error: Date of extraction must be in format %YYYY-%mm-%dd")
        sys.exit(1)
    
    if args.output_file is None:
        args.output_file = Path(f'data/{"somatic" if args.somatic else "germline"}_{args.variant_type}_uploaded_annotated_{date_of_interest}.tsv')
    else: #check it is saved in 'data' folder
        if args.output_file.split("/")[0] != "data":
            args.output_file = f'data/{args.output_file}'
    annotate_file(args.summary_files, args.input_file, args.reference_file, args.output_file, args.variant_type, date_of_interest)
