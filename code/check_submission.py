import argparse
import requests
from pathlib import Path
import helper_functions

def parse_arguments():
    parser = argparse.ArgumentParser(description="Specify Submission ID of the batch to be check whether the upload has been succesfull and return the summary json.")
    parser.add_argument(
        "-i",
        "--submission_id",
        help="Submission ID of batch of variants to be checked (SUBnnnnnnnn).",
        required=True,
    )
    parser.add_argument(
        "-k",
        "--key",
        help='Path to "clinvar.key" file',
        type=Path,
        default=Path("./clinvar.key"),
    )
    parser.add_argument(
        "-n",
        "--dryrun",
        help="If batch was submitted to the dry run URL",
        action="store_true",
    )
    return parser

if __name__ == "__main__":
    parser = parse_arguments()
    args = parser.parse_args()
    if args.dryrun:
        check_url = f"{helper_functions.TEST_URL}/{args.submission_id}/actions/"
    else:
        check_url = f"{helper_functions.SUB_URL}/{args.submission_id}/actions/"

    api_key = args.key.read_text()
    headers = {
        "Content-Type": "application/json",
        "SP-API-KEY": api_key
        }

    sub_response = requests.get(
        url = check_url,
        headers = headers
    )

    print(sub_response.json())
