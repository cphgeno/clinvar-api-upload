import copy
import csv
from datetime import date
import os
from re import sub

BASE_PAIRS = {"A": "T", "C": "G", "G": "C", "T": "A"}
SUB_URL = "https://submit.ncbi.nlm.nih.gov/api/v1/submissions"
TEST_URL = "https://submit.ncbi.nlm.nih.gov/apitest/v1/submissions"

SOMATIC_CLASSIFICATION_MAPPING = {"Pathogenic": "Oncogenic", "Likely pathogenic": "Likely oncogenic", "Uncertain significance": "Uncertain significance", "Likely benign":"Likely benign", "Benign": "Benign"}

def get_reverse_strand(sequence):
    """Function to obtain the reverse-complement of a genomic sequence.

    Parameters:
        sequence: template sequence from which to obtain the new strand.

    Return:
        Reverse-complement of the sequence"""

    comp = [BASE_PAIRS[nt] for nt in list(sequence)]
    return "".join(comp[::-1])


def fix_coordinates(sample):
    """Function to modify the variant's dictionary when \"hgvs c.\" is not provided.

    In the case of a SNP, the start coordinate must be increased by one to match the stop coordinate.

    In the case of an insertion, the stop coordinate must be increased by one.

    In the case of a deletion, the start coordinate must be increased by one to match the reference sequence.

    Parameters:
        sample: sample to be modified.

    Return:
        Modified sample"""

    ref_allele = sample["variantSet"]["variant"][0]["chromosomeCoordinates"][
        "referenceAllele"
    ]
    alt_allele = sample["variantSet"]["variant"][0]["chromosomeCoordinates"][
        "alternateAllele"
    ]
    if (
        sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["start"]
        == sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["stop"]
    ) or (
        int(sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["stop"])
        - int(sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["start"]) == 1
    ):
        if ref_allele == "-":
            sample["variantSet"]["variant"][0]["variantType"] = "Insertion"
            sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["stop"] += 1
    else:
        if alt_allele == "-":
            sample["variantSet"]["variant"][0]["variantType"] = "Deletion"
            sample["variantSet"]["variant"][0]["chromosomeCoordinates"]["start"] += 1
        else:
            if len(ref_allele) == len(alt_allele) == 1:
                sample["variantSet"]["variant"][0]["chromosomeCoordinates"][
                    "start"
                ] += 1
    return sample


def coordinates_check(variants_list):
    """Function used to make sure data is in the correct format for upload.
    If variant has hgvs c. value then no need to check as already in ocrrect format for the previous cleaning step.
    
    Parameters:
        variants_list: variants to be checked before upload.

    Return:
        Modified variants list."""
    final_list = []
    for sample in variants_list:
        if "chromosomeCoordinates" in sample["variantSet"]["variant"][0].keys():
            final_list.append(fix_coordinates(sample))
        else:
            final_list.append(
                sample
            )
    return final_list


def fix_hgvs(hgvs, alt):
    """Function to modify the HGVS syntax when uncertain.

    If parenthesis are present in the HGVS, this syntax is not accepted by the API, hence the number in parenthesis, representing the length of the sequence, is replaced by the sequence itself.

    Parameters:
        hgvs: HGVS of the variant of interest.
        alt: alternate allele of the variant of interest.

    Return:
        Fixed HGVS."""

    modified_hgvs = []
    if "(" in hgvs:
        values = hgvs.split(",")
        if len(values) == 2:
            value1 = sub("[\(].*?[\)]", get_reverse_strand(alt), values[0])
            modified_hgvs.append(value1)
            value2 = sub("[\(].*?[\)]", alt, values[1])
            modified_hgvs.append(value2)
        else:
            modified_hgvs.append(sub("[\(].*?[\)]", get_reverse_strand(alt), values[0]))
        #print(f'Old hgvs c.: {hgvs} \n New hgvs c.: {",".join(modified_hgvs)}')
    return ",".join(modified_hgvs)


def check_hgvs(sample):
    """Function to check whether the HGVS is present and its syntax correct.

    Parameters:
        sample: sample to be checked.

    Return:
        Boolean value: "True" if hgvs is present, "False" if not."""

    if sample["hgvs c."]:
        if "(" in sample["hgvs c."]:
            sample["hgvs c."] = fix_hgvs(
                sample["hgvs c."], sample["Ref/Alt"].split("/")[1]
            )
        return True
    return False


def extract_variants(hgvs_list, variant):
    """Function to extract the variants when more than one HGVS is present.

    Each HGVS provided is separated together with the gene associated to it and they are placed into a new variant.
    All other data is copied as it does not change.

    Parameters:
        hgvs_list: list of the multiple HGVS found in the variant.
        variant: variant of interest.

    Return:
        List of single extracted variants from the inputted one."""

    variants_list = []
    genes_list = variant["variantSet"]["variant"][0]["gene"][0]["symbol"].split(",")
    for i in range(len(hgvs_list)):
        if not hgvs_list[i].startswith("NM"):
            continue
        else:
            single_variant = copy.deepcopy(variant)
            single_variant["variantSet"]["variant"][0]["hgvs"] = hgvs_list[i]
            single_variant["variantSet"]["variant"][0]["gene"][0][
                "symbol"
            ] = genes_list[i]
            variants_list.append(single_variant)
    return variants_list


def multiple_hgvs(sample):
    """Function to identify variants possessing more than one HGVS.

    When a variant presents multiple HGVS values, these are separated into new variants.
    Usually the case when two genes, on opposite strands, are affected by the variant.
    Data is untouched if only a single HGVS value is present.

    Parameters:
        sample: variant of interest to be checked.

    Return:
        The separated variants if multiple HGVS present or "None" otherwise."""

    hgvs_list = sample["variantSet"]["variant"][0]["hgvs"].split(",")
    if len(hgvs_list) == 1:
        return None
    else:
        single_variants = extract_variants(hgvs_list, sample)
        return single_variants


def convert_variant(variants_entries, date_of_upload, somatic_flag, to_update=False):
    """Function to convert the variant into the form required for ClinVar upload.

    The initial dictionary read from the input file is converted to another dictionary reflecting the keys specified by the Submission API (https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/).

    Parameters:
        variants_entries: variants to be formatted according to ClinVar API schema.
        date_of_upload: date of upload to be included in submission of haplotype.
        to_update: boolean value to determine whether variants are to be updated or just uploaded.

    Returns:
        variants_list: list of variants correctly formatted."""

    variants_list = []
    for variant in variants_entries:
        variant["Classification"] = variant["Classification"].capitalize()
        if variant["Classification"] == "Unknown significance":
            variant["Classification"] = "Uncertain significance"  # ClinVar requirement
        if somatic_flag:
            variant["Classification"] = SOMATIC_CLASSIFICATION_MAPPING[variant["Classification"]]
        
        condition_second_field = (
            {"id": "C0027651"} if somatic_flag else
            {"name": (
                "not specified"
                if variant["Classification"] in ["Benign", "Likely benign", "Uncertain significance"]
                else "not provided"
            )}
        )
        
        sample_variant = {
            ("germlineClassification" 
             if not somatic_flag
             else "oncogenicityClassification"): {
                ("germlineClassificationDescription"
                 if not somatic_flag
                 else "oncogenicityClassificationDescription"): variant["Classification"],
                "dateLastEvaluated": date_of_upload,
            },
            "conditionSet": {
                "condition": [
                    {
                        "db": "MedGen",
                        **condition_second_field
                    }
                ]
            },
            "observedIn": [
                {
                    "affectedStatus": (
                        "unknown" if not somatic_flag # we choose 'unknown' as default
                        else "yes" 
                    ),
                    "alleleOrigin": (
                        "germline" if not somatic_flag
                        else "somatic"
                    ),
                    "collectionMethod": "clinical testing",  # we choose 'clinical testing' as default
                }
            ],
            "recordStatus": "novel" if not to_update else "update",
            "variantSet": {
                "variant": [
                    {
                        "gene": [{"symbol": variant["Gene Names"]}],
                        (
                            "hgvs"
                            if check_hgvs(variant)
                            else "chromosomeCoordinates"
                        ): (
                            variant["hgvs c."]
                            if variant["hgvs c."]
                            else {
                                "assembly": "hg38",
                                "chromosome": variant["#Chromosome"],
                                "start": int(variant["Start"]),
                                "stop": int(variant["Stop"]),
                            }
                        ),
                    }
                ]
            },
        }
        if to_update:  # if variants to be updated then SCV field required
            sample_variant["clinvarAccession"] = variant["SCV"]
        # reference and alternate alleles are only to be provided if variant is smaller than 50nt
        if (
            "chromosomeCoordinates" in sample_variant["variantSet"]["variant"][0].keys()
            and int(variant["Stop"]) - int(variant["Start"]) <= 50
        ):
            sample_variant["variantSet"]["variant"][0]["chromosomeCoordinates"][
                "referenceAllele"
            ] = variant["Ref/Alt"].split("/")[0]
            sample_variant["variantSet"]["variant"][0]["chromosomeCoordinates"][
                "alternateAllele"
            ] = variant["Ref/Alt"].split("/")[1]
        variants_list.append(sample_variant)
    return variants_list


def format_haplo_variants(variants_list):
    """Function to obtain the list of variants that are part of the haplotype and need to be included in the haplotypes submission.
    The syntax on VarSeq should result in the hgvs c. to be present for all of them.
    
    Parameters:
        variants_list: list of variants extracted from the VarSeq haplotype's syntax.
    
    Returns:
        output_variants: formatted list of variants for haplotype submission.
    """

    output_variants = []
    for hgvsc in variants_list.split(", "):
        formatted_variant = {}
        formatted_variant["hgvs"] = hgvsc
        output_variants.append(formatted_variant)
    return output_variants


def convert_haplotype(haplo_entries, date_of_upload, somatic_flag, to_update):
    """Function to convert the haplotypes into the form required for ClinVar upload.

    The initial dictionary read from the input file is converted to another dictionary reflecting the keys specified by the Submission API (https://www.ncbi.nlm.nih.gov/clinvar/docs/api_http/).

    Parameters:
        haplo_entries: haplotypes to be formatted according to ClinVar API schema.
        date_of_upload: date of upload to be included in submission of haplotype.

    Returns:
        haplotypes_list: list of haplotypes correctly formatted for submission."""

    haplotypes_list = []

    for haplo in haplo_entries:
        if to_update:
            associated_variants, haplo_hgvsc, haplo_hgvsp, haplo_classification, upload_type, last_edited_date, haplo_SCV = haplo_entries[haplo].strip(" ").split("; ")
        else:
            associated_variants, haplo_hgvsc, haplo_hgvsp, haplo_classification, upload_type, last_edited_date = haplo_entries[haplo].strip(" ").split("; ")
        haplo_classification = haplo_classification.capitalize()
        if haplo_classification == "Unknown significance":
            haplo_classification = "Uncertain significance"
        if somatic_flag:
            haplo_classification = SOMATIC_CLASSIFICATION_MAPPING[haplo_classification]
        
        condition_second_field = (
            {"id": "C0027651"} if somatic_flag else
            {"name": (
                "not specified"
                if haplo_classification in ["Benign", "Likely benign", "Uncertain significance"]
                else "not provided"
            )}
        )

        sample_variant = {
            ("germlineClassification" 
             if not somatic_flag
             else "oncogenicityClassification"): {
                ("germlineClassificationDescription"
                 if not somatic_flag
                 else "oncogenicityClassificationDescription"): haplo_classification,
                "dateLastEvaluated": date_of_upload,
            },
            "conditionSet": {
                "condition": [
                    {
                        "db": "MedGen",
                        **condition_second_field
                    }
                ]
            },
            "observedIn": [
                {
                    "affectedStatus": (
                        "unknown" if not somatic_flag
                        else "yes" 
                    ),
                    "alleleOrigin": (
                        "germline" if not somatic_flag
                        else "somatic"
                    ),
                    "collectionMethod": "clinical testing",
                }
            ],
            "recordStatus":"novel" if not to_update else "update",
            "haplotypeSet": {
                "hgvs": haplo_hgvsc,
                "variants": format_haplo_variants(associated_variants)
            },
        }
        if to_update:  # if variants to be updated then SCV field required
            sample_variant["clinvarAccession"] = haplo_SCV
        haplotypes_list.append(sample_variant)
    return haplotypes_list

def write_haplotypes_file(haplotypes_dict, file_path, date_of_upload):
    """Function for writing the list of uploaded haplotypes to a txt file.
    This file is created if not already present, otherwise the newly uploaded haplotypes are appended to the existing file.

    Parameters:
        haplotypes_dict: dictionary of haplotypoes to be stored in the file.
        file_path: output's filepath.
        date_of_upload: date of upload needed for annotation.
    """
    lines = []
    for haplotype in haplotypes_dict:
        associated_variants, haplo_hgvsc, haplo_hgvsp, haplo_classification, upload_type, last_edited = haplotypes_dict[haplotype].strip(" ").split("; ")[:6] #if to be updated skip SCV and add later again with annotation
        lines.append([haplo_hgvsc, haplo_classification, associated_variants, haplo_hgvsp, date_of_upload])
    if os.path.exists(file_path):
        with open(file_path, 'a', newline='', encoding='utf-8') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerows(lines)
    else:
        with open(file_path, 'w', newline='', encoding='utf-8') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            print("File created")
            header = ["hgvs c.", "Classification", "Variants", "hgvs p.", "Last Edited"]
            writer.writerow(header)
            writer.writerows(lines)
    return

def format_deletion(scv_entries):
    """Function responsible for formatting the variants that need to be deleted by the SCV accession number.
    
    Parameters:
        scv_entries: list of SCV accession numbers of the variants to be deleted.
    
    Output:
        variants_list: list of formatted SCV accession numebrs for submission to the API.
    """
    variants_list = {
        "accessionSet":
                [{"accession": scv_accession} for scv_accession in scv_entries]
    }
    return variants_list
