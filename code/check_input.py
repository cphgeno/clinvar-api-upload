import csv, copy
import re
from collections import defaultdict


def separate_variants(two_variants):
    """Function to separate the lines with three alleles into two separate variants.

    Common in repetitive regions that two alternate alleles are reported for a single variant in the same line.
    ...

    Parameters:
        two_variants: dicitonary element containing the variant to be split.

    Return:
        the two single variants."""

    v1 = copy.deepcopy(two_variants)
    v1["Ref/Alt"] = "/".join(two_variants["Ref/Alt"].split("/")[:2])
    v2 = copy.deepcopy(two_variants)
    v2["Ref/Alt"] = "/".join(two_variants["Ref/Alt"].split("/")[::2])
    if "," in two_variants["hgvs c."]:
        v1["hgvs c."], v2["hgvs c."] = two_variants["hgvs c."].split(",")
    else:
        v2[
            "hgvs c."
        ] = ""  # if only one present it describes only first of two variants
    return v1, v2


def correct_format(ref, alt, start, stop, cursor, step):
    """Function to reduce the variants to the minimal writing of alleles and consequently modifying the genomic position.

    It compares the strings of reference and alternate alleles until one of them ends or there is a difference.
    The direction of comparison is determmined by the cursor and step values.

    Parameters:
        ref: reference allele.
        alt: alternate allele.
        start: start position of variant.
        stop: stop position of variant.
        cursor: for indexing of the alleles (strings), either 0 or -1.
        step: indicating which direction the comparison of the strings will happen.

    Return:
        New parameters for reference, alternate alleles, start and stop positions."""

    subset = False
    while abs(cursor) < len(ref) and abs(cursor) < len(alt):
        if ref[cursor] != alt[cursor]:
            subset = (
                True  # if only a subset of both alleles is in common, iteration stopped
            )
            break
        cursor += step
    if subset and cursor < -1:  # suffix in common
        # cursor has to be increased by one for proper indexing and location
        ref = ref[: cursor + 1]  # extract new reference allele
        alt = alt[: cursor + 1]  # extract new alternate allele
        stop = str(
            int(stop) + cursor + 1
        )  # stop position is updated base on length of suffix
        if ref[0] == alt[0]:  # if also prefix in common the function is recalled
            return correct_format(ref, alt, start, stop, cursor=0, step=1)
    elif subset and cursor > 0:  # prefix in common
        ref = ref[cursor:]  # update alleles
        alt = alt[cursor:]
        start = str(
            int(start) + cursor
        )  # start position updated based on length of prefix
    # if whole reference is subset of alternate allele
    elif abs(cursor) > len(ref):  # insertion
        if cursor > 0:
            alt = alt[cursor:]
            start = str(int(start) + cursor)
        else:
            # N.B. again cursor to be increased by 1 when negative
            alt = alt[: cursor + 1]
            stop = str(int(stop) + cursor + 1)
        ref = ["-"]  # reference allele becomes empty as is insertion
    # if whole alternate is subset of reference allele
    # same reasoning as for the insertion case above
    elif abs(cursor) > len(alt):  # deletion
        if cursor > 0:
            ref = ref[cursor:]
            start = str(int(start) + cursor)
        else:
            ref = ref[: cursor + 1]
            stop = str(int(stop) + cursor + 1)
        alt = ["-"]
    # when cursor is same as allele length
    # extra cases where while loop ended without finding subset
    elif abs(cursor) == len(ref):
        if not ref[cursor] == alt[cursor] and cursor < 0:
            # if first nt of alleles are different
            ref = ref[0]
            if (
                cursor == -1
            ):  # if while loop is not even executed (ref allele has length 1)
                alt = alt  # no change to the reference allele
            else:
                alt = alt[: cursor + 1]
            stop = str(int(stop) + cursor + 1)
        elif cursor < 0:  # if they are the same, i.e. insertion
            if (
                cursor == -1
            ):  # if while loop is not even executed (ref allele has length 1)
                alt = alt[:cursor]  # no need to increase cursor in this case
            else:
                alt = alt[:cursor]
                stop = str(int(stop) + cursor)  # no need to increase cursor
            ref = ["-"]
        else:  # if cursor >0 then it is an insertion
            # comparison of alleles has already been done because indexing starts at 0, not 1
            alt = alt[cursor:]
            ref = ["-"]
            start = str(int(start) + cursor)
    # same reasoning as above
    elif abs(cursor) == len(alt):
        if not ref[cursor] == alt[cursor]:
            if cursor > 0:
                ref = ref[cursor:]
                alt = alt[-1]
                start = str(int(start) + cursor)
            else:
                ref = ref[: cursor + 1]
                alt = alt[0]
                stop = str(int(stop) + cursor)
        else:
            if cursor > 0:
                ref = ref[cursor:]
                start = str(int(start) + cursor)
            else:
                ref = ref[:cursor]
                stop = str(int(stop) + cursor)
            alt = ["-"]
    return "/".join(["".join(ref), "".join(alt)]), start, stop


def check_variant(variant, reference_file):
    """Function to check whether the variant is already present in the file after having obtained the alleles in the correct format.

    Parameters:
        variant: variant to be checked.

    Return:
        Boolean value: "True" if the variant is not present in the file and to be added to it, "False" oterwise."""

    # modify alleles if there is overlap
    ref, alt = variant["Ref/Alt"].split("/")
    if (
        not (alt != "-" and ref == "-")
        and not (ref != "-" and alt == "-")
        and not (len(ref) == len(alt) == 1)
    ):
        if (
            ref[-1] == alt[-1]
        ):  # check for common suffix, and also for possible prefix afterwards
            variant["Ref/Alt"], variant["Start"], variant["Stop"] = correct_format(
                ref, alt, variant["Start"], variant["Stop"], cursor=-1, step=-1
            )
            variant[
                "hgvs c."
            ] = ""  # delete hgvs c. as will be wrong, possibly to be generated manually
        elif ref[0] == alt[0]:  # check for common prefix
            variant["Ref/Alt"], variant["Start"], variant["Stop"] = correct_format(
                ref, alt, variant["Start"], variant["Stop"], cursor=0, step=1
            )
            variant[
                "hgvs c."
            ] = ""  # delete hgvs c. as will be wrong, possibly to be generated manually
    variant_check = "\t".join(
        [variant["#Chromosome"], variant["Start"], variant["Stop"], variant["Ref/Alt"]]
    )

    # check if variant already present in file being written
    with open(reference_file, encoding="utf-8", mode="r+") as check_outfile:
        outfile_reader = csv.DictReader(check_outfile, delimiter="\t")
        for element in outfile_reader:
            if variant["hgvs c."]:
                if variant["hgvs c."] == element["hgvs c."]:
                    return False
            check_line = "\t".join(
                [
                    element["#Chromosome"],
                    element["Start"],
                    element["Stop"],
                    element["Ref/Alt"],
                ]
            )
            if variant_check == check_line:
                return False
    return True


def extract_variants(variant):
    """Function to extract the variants when more than one HGVS is present.

    Because of the format ClinVar requires, only one hgvs c. is kept and the
    others arediscarded as it associates the variant automatically to the other
    genes it affects (both strands).

    Parameters:
        hgvs_list: list of the multiple HGVS found in the variant.
        variant: variant of interest.

    Return:
        First variant from inputted list that contains both correct hgvs c. syntax and gene name."""

    hgvs_list = variant["hgvs c."].split(",")
    genes_list = variant["Gene Names"].split(",")
    for i in range(len(hgvs_list)):
        if not hgvs_list[i].startswith("NM"):
            continue
        else:
            single_variant = copy.deepcopy(variant)
            single_variant["hgvs c."] = hgvs_list[i]
            single_variant["Gene Names"] = genes_list[i]
            return single_variant


def remove_duplicate_rows(file_path, delimiter='\t'):
    """
    Function to remove potential duplicated rows in the inputted file.

    Parameters:
        file_path: initial tsv file containing the variants.
    
    Return:
        unique_rows: list of unique dictionary variants to be parsed through.
        fieldnames: header from the input file to be used for cleaned file writing.
    """

    unique_rows = []
    seen = defaultdict(int)  # Dictionary to track occurrences
    
    with open(file_path, 'r', newline='', encoding='utf-8') as file:
        csv_reader = csv.DictReader(file, delimiter=delimiter)
        fieldnames = csv_reader.fieldnames
        rows = list(csv_reader)
        
        for row in rows:
            key = (row['hgvs c.'], '<->'.join([row['#Chromosome'], row['Start'], row['Stop'], row['Ref/Alt']]))
            seen[key] += 1
            if seen[key] > 1:
                continue
            else:
                unique_rows.append(row)
    
    print(f'Number or duplicated variants removed from input file: {len(rows) - len(unique_rows)}')
    
    return unique_rows, fieldnames


def clean_data(input_file):
    """Function that allows to clean the inputted tsv file.
    The function parses through the variants and executes three main functionalities:
        - gets rid of variants that report more than one alternate allele;
        - removes duplicate variants;
        - identifies variants for which more than one hgvs c. is present and cleans them;
        - extracts the haplotypes information and removes the individual variants part of it if not to be uploaded individually.

    Parameters:
        input_file: initial tsv file containing the variants.

    Return:
        output_file: modified/cleaned tsv file with the variants.
        output_file_haplo: tsv file with haplotypes only information for parsing to later functions."""

    file_basename = input_file.split("/")[-1].split(".")[0]
    output_file = f'temp/{file_basename}_cleaned.txt'
    output_file_haplo = f'temp/{file_basename}_haplotypes_cleaned.txt'

    variants_reader, fieldnames = remove_duplicate_rows(input_file)

    with open(output_file, encoding="utf-8", mode="w", newline="") as outfile:
        csv_writer = csv.DictWriter(
            outfile, fieldnames = fieldnames, delimiter="\t"
        )
        csv_writer.writeheader()

        triple_variants = []
        multiple_hgvs = []
        haplotypes_dict = {}
        for variant in variants_reader:
            if len(variant["Ref/Alt"].split("/")) == 3:
                triple_variants.append(variant)
            elif "," in variant["hgvs c."]:
                if not "[MERGE:" in variant["Notes"]:
                    multiple_hgvs.append(variant)
            elif "[MERGE:" in variant["Notes"]:
                merge_info = re.search(r'\[MERGE:(.+)\]', variant["Notes"]).group(1)
                merge_info_mod = merge_info.replace("&gt;", ">")
                associated_variants, haplo_hgvsc, haplo_hgvsp, haplo_classification, upload_type = merge_info_mod.strip(" ").split("; ")
                last_edited_date = variant["Last Edited"]
                if haplo_hgvsc not in haplotypes_dict:
                    # haplotypes_dict[haplo_hgvsc] = "; ".join([merge_info_mod, last_edited_date])
                    haplotypes_dict[haplo_hgvsc] = variant # does not account for difference in last edited for both variants, possibly to be changed later
                if upload_type == "individual-merged":
                    # variant to be uploaded individually AND as part of haplotype
                    csv_writer.writerow(variant)
            else:
                csv_writer.writerow(variant)
        new_variants = []
    
    if multiple_hgvs:
        for variant in multiple_hgvs:
            new_variant = extract_variants(variant)
            new_variants.append(new_variant)

        with open(
            output_file, encoding="utf-8", mode="a", newline=""
        ) as check_outfile:
            csv_writer = csv.DictWriter(
                check_outfile, fieldnames= fieldnames, delimiter="\t"
            )
            csv_writer.writerows(new_variants)
    
    if triple_variants:
        for variant in triple_variants:
            v1, v2 = separate_variants(variant)
            if check_variant(v1, output_file):
                with open(
                    output_file, encoding="utf-8", mode="a", newline=""
                ) as check_outfile:
                    csv_writer = csv.DictWriter(
                        check_outfile, fieldnames= fieldnames, delimiter="\t"
                    )
                    csv_writer.writerow(v1)
            if check_variant(v2, output_file):
                with open(
                    output_file, encoding="utf-8", mode="a", newline=""
                ) as check_outfile:
                    csv_writer = csv.DictWriter(
                        check_outfile, fieldnames= fieldnames, delimiter="\t"
                    )
                    csv_writer.writerow(v2)

    # write haplotypes file, same format as variants for upload with same function
    if not haplotypes_dict:
        output_file_haplo = None
    else:
        with open(output_file_haplo, encoding="utf-8", mode="w") as haplo_tsv:
            csv_writer = csv.DictWriter(
                haplo_tsv, fieldnames=['hgvs c.', 'Classification', 'Notes', 'hgvs p.', 'Gene Names', 'Kriterier', 'Last Edited'], delimiter="\t"
            )
            csv_writer.writeheader()

            for haplo_hgvs, variant in haplotypes_dict.items():
                merge_info = re.search(r'\[MERGE:(.+)\]', variant["Notes"]).group(1)
                merge_info_mod = merge_info.replace("&gt;", ">")
                associated_variants, haplo_hgvsc, haplo_hgvsp, haplo_classification, upload_type = merge_info_mod.strip(" ").split("; ")

                csv_writer.writerow({
                    'hgvs c.': haplo_hgvs, 
                    'Classification': haplo_classification,
                    'Notes': variant['Notes'],
                    'hgvs p.': haplo_hgvsp,
                    'Gene Names': variant['Gene Names'],
                    'Last Edited': variant['Last Edited']
                })

    return output_file, output_file_haplo

