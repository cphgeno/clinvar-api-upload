import argparse
import csv
import json
import os
import sys
import requests
from pathlib import Path
import helper_functions, check_input, check_for_updates
from datetime import datetime

os.makedirs(".\\temp\\")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Provide data to be uploaded and reference file.")
    parser.add_argument(
        "-f",
        "--file",
        help="Path to VarSeq Assessment Catalog exported tsv",
        required=True
    )
    parser.add_argument(
        "-k",
        "--key",
        help='Path to "clinvar.key" file',
        type=Path,
        default=Path("./clinvar.key")
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        help="To run data validation testing (dry run)",
        action="store_true"
    )
    parser.add_argument(
        "-s",
        "--status",
        help='Set to "update" if the submitted variants are already present on ClinVar and need to be updated (default = "update")',
        choices=["novel", "update", "delete"],
        default="update"
    )
    parser.add_argument(
        "--reference_variants",
        help="Reference file, containing variants that have already been uploaded, required for update, can be found on the GitLab repository",
        required=True
    )
    parser.add_argument(
        "--reference_haplotypes",
        help="Reference file, containing haplotypes that have already been uploaded, required for update, can be found on the GitLab repository",
        required=True
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        help="Specify maximum batch size for submisison (default = 500)",
        type=int,
        default=500
    )
    parser.add_argument(
        "--date_of_extraction",
        help="Date of extraction of variants from VarSeq (API upload date does not reflect this as always done with some delay) - important for data anonimisation",
        required=True
    )
    return parser


def submit_request(list_of_variants, submission_url, headers, deletion = False):
    """Function to post the request for submission of each batch of variants.
    If the deletion argument is set to True, the function posts the request for the deletion of the submitted variants.

    Parameters:
        list_of_variants (list): list of the batch of variants to be submitted.
        submission_url: submission url for dry or live run (based on -n argument).
        headers: headers with extra info needed for the submission.
        deletion: Boolean value for submission of deletion of variants.

    Return:
        Response object (dict) from the attempted submission."""

    if not deletion:
        content = {
            "clinvarSubmission": list_of_variants,
            "clinvarSubmissionReleaseStatus": "public",
            "assertionCriteria":{
                "db": "PubMed",
                "id": "" # INSERT: PubMed id of the article for assertion criteria
            }
        }
    else:
        content = {"clinvarDeletion": list_of_variants}
    data = {
        "actions": [
            {"type": "AddData", "targetDb": "clinvar", "data": {"content": content}}
        ]
    }
    response = requests.post(url=submission_url, headers=headers, data=json.dumps(data))
    return response


def submit_variants(variants_list, variant_status, haplo = False):
    """Wrapper function responible for uploading the novel or updated variants.

    Parameters:
        varinats_list: list of variants to be uploaded.
        variant_status: string indicating whether the variants are novel or to be updated.
        haplo: Boolean value to define whether variants or haplotypes are to be uploaded."""

    status_choices = ["novel", "update", "delete"]
    if variant_status not in status_choices:
        raise ValueError("Invalid sim type. Expected one of: %s" % status_choices)
    if args.dryrun:
        submission_url = helper_functions.TEST_URL
    else:
        submission_url = helper_functions.SUB_URL
    
    api_key = args.key.read_text()
    headers = {"Content-Type": "application/json", "SP-API-KEY": api_key}

    if variant_status == "delete":
        response = submit_request(variants_list, submission_url, headers, deletion = True)
        print(f"Response ({variant_status}): {response.status_code}")
        sub_id = response.text.strip("{}").split(":")[1].strip('"')
        print(f"Submission id: {sub_id}")
        print(response.headers)
        return
    
    print(f"[INFO]: Submitting variants for upload ({variant_status})...")
    responses = []
    variants_number = len(variants_list)
    variants_list = [
        variants_list[i : i + args.batch_size]
        for i in range(0, len(variants_list), args.batch_size)
    ]
    for batch in variants_list:
        responses.append(submit_request(batch, submission_url, headers))
    # -------on-screen message
    print(f'{"Live run" if args.dryrun == False else "Dry run"} at: {submission_url}')
    print(f"Number variants submitted: {variants_number}")
    summaries_list = []

    for count, response in enumerate(responses, 1):
        print(f"Response ({variant_status}): {response.status_code}")
        sub_id = response.text.strip("{}").split(":")[1].strip('"')
        print(f"Submission id: {sub_id}")
        print(response.text)

        if response.status_code > 201:
            with open(
                f"Submission_errors_{variant_status}_{count}.txt", mode="w"
            ) as output_file:
                output_file.write(str(response.headers))
                output_file.write("\n")
                output_file.write(response.text)
        else:
            data = json.dumps(variants_list, indent=2)
            with open(f'temp\\{sub_id}_{"variants" if not haplo else "haplotypes"}_{variant_status}.json', mode="w") as output_file:
                header = dict(response.headers)
                header["Total number of variants submitted"] = variants_number
                output_file.write(json.dumps(header, indent=2))
                output_file.write(data)
            report_json = f'temp\\{sub_id}-summary-report.json {"variants" if not haplo else "haplotypes"}'
        if haplo and variant_status=="novel":
            helper_functions.write_haplotypes_file(haplos_novel, "temp\haplotypes_uploaded.txt", date_of_extraction)
        summaries_list.append(report_json)
        file_path = f'temp\summaries_list_{variant_status}.txt'
        with open(file_path, mode = 'a' if os.path.exists(file_path) else 'w') as summaries_file:
            #writes file with summary jsons path to be used for annotation
            for json_file in summaries_list:
                summaries_file.write(json_file)
                summaries_file.write("\n")
        print("[INFO]: Variant submission to ClinVar API completed!")
        print("------IMPORTANT! Remember to annotate the data after successful upload------ \n \n")

# ----------main execution-------------------
if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()

    if args.status == "update":
        if args.reference_variants is None:
            parser.error("The --reference_variants flag is required to update the variants")
        if args.reference_haplotypes is None:
            parser.error("The --reference_haplotypes flag is required to update the haplotypes")
    if len(args.date_of_extraction) == 10:
        try:
            datetime.strptime(args.date_of_extraction, '%Y-%m-%d')
            # date_check = bool(date_parser.parse(args.date_of_extraction))
            date_of_extraction = args.date_of_extraction
        except ValueError:
            print("Error: Date of extraction must be in format %YYYY-%mm-%dd")
            sys.exit(1)
    else:
        print("Error: Date of extraction must be in format %YYYY-%mm-%dd")
        sys.exit(1)
    
    if args.status == "delete":
        print("[INFO]: Submitting variants to be deleted...")
        with open(args.file, encoding="utf-8", mode="r") as input_tsv:
            csv_reader = csv.DictReader(input_tsv, delimiter="\t")
            scv_list = [variant["SCV"] for variant in csv_reader]
        variants_list = helper_functions.format_deletion(scv_list)
        submit_variants(variants_list, variant_status = "delete")
    else: #in case of upload/update
        print(f"[INFO]: Cleaning data from {args.file}...")
        input_file, haplotypes_dict = check_input.clean_data(args.file)
        if args.status == "update":
            variants_novel, variants_update = check_for_updates.compare_variants(
                input_file, args.reference_variants, date_of_extraction
            )

            print(f"[INFO]: Parsing through the novel variants...")
            if len(variants_novel) != 0:
                print(len(variants_novel))
                # print(variants_novel)
                variants_novel_formatted = helper_functions.convert_variant(variants_novel, date_of_extraction)
                variants_novel_formatted = helper_functions.coordinates_check(variants_novel_formatted)
                submit_variants(variants_novel_formatted, variant_status="novel")
            else:
                print("[INFO]: No novel variants to be uploaded.")
            
            print(f"[INFO]: Parsing through the variants to be updated...")
            if len(variants_update) != 0:
                print(len(variants_update))
                # print(variants_update)
                variants_update_formatted = helper_functions.convert_variant(variants_update, date_of_extraction, to_update=True)
                variants_update_formatted = helper_functions.coordinates_check(variants_update_formatted)
                submit_variants(variants_update_formatted, variant_status="update")
            else:
                print("[INFO]: No variants to be updated.")
            
            print(f"[INFO]: Looking at haplotypes...")
            if haplotypes_dict:
                haplos_novel, haplos_update = check_for_updates.compare_haplotypes(haplotypes_dict, args.reference_haplotypes)
                if len(haplos_novel) != 0:
                    haplotypes_upload_formatted = helper_functions.convert_haplotype(haplos_novel, date_of_extraction)
                    submit_variants(haplotypes_upload_formatted, variant_status="novel", haplo = True)
                else:
                    print("[INFO]: No novel haplotypes to be uploaded")
                if len(haplos_update) != 0:
                    haplotypes_update_formatted = helper_functions.convert_haplotype(haplos_update, date_of_extraction)
                    submit_variants(haplotypes_update_formatted, variant_status="update", haplo = True)
                else:
                    print("[INFO]: No haplotypes to be updated.")
        else:
            print(f"[INFO]: Parsing through the novel variants...")
            with open(input_file, encoding="utf-8", mode="r") as tsv:
                csv_reader = csv.DictReader(tsv, delimiter="\t")
                input_entries = []
                for variant in csv_reader:
                    input_entries.append(variant)
            variants_list = helper_functions.convert_variant(input_entries, date_of_extraction)
            submit_variants(variants_list, variant_status="novel")
