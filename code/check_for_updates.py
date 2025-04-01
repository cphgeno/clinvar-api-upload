import csv

def read_file(file_variants):
    """Function to read the variants into a dictionary.

    Parameteres:
        file_variants: filepath to latest tsv file containing the variants.

    Returns:
        file_entries: dictionary containing all variants from the corresponding file."""

    with open(file_variants, encoding="utf-8", mode="r") as f:
        file_reader = csv.DictReader(f, delimiter="\t")
        file_entries = []
        for new_variant in file_reader:
            file_entries.append(new_variant)
    return file_entries

def compare_variants(new_variants, old_variants, date_of_upload, haplo = False):
    """Function to compare the new variants from the input file to the reference file.

    If a variant in the input file is present in the reference file the "Last Edited" date is checked.
    If the date is different, the variants needs to be updated, if the same the latest version of it is already present on ClinVar.
    If the variant is not present in the reference file it means it is novel and has t be uploaded from scratch.

    Parameters:
        new_variants: list containing variants from the input file.
        old_variants: list containing variants from the reference file.
        date_of_upload: date of upload to modify in the variant's field, if uploaded, for future comparison.

    Return:
        variants_upload: list containing variants to be uploaded as novel.
        variants_update: list containing variants to be updated."""

    new_entries = read_file(new_variants)
    old_entries = read_file(old_variants)
    variants_update = []
    variants_upload = []
    for new_variant in new_entries:
        if not haplo:
            new_variant_coords = "\t".join(
                [
                    new_variant["#Chromosome"],
                    new_variant["Start"],
                    new_variant["Stop"],
                    new_variant["Ref/Alt"],
                ]
            )
        to_update = False
        up_to_date = False
        for old_variant in old_entries:
            if not haplo:
                old_variant_coords = "\t".join(
                    [
                        old_variant["#Chromosome"],
                        old_variant["Start"],
                        old_variant["Stop"],
                        old_variant["Ref/Alt"],
                    ]
                )
            if (
                new_variant["hgvs c."] == old_variant["hgvs c."] != "") or (
                not haplo and new_variant_coords == old_variant_coords
            ):
                if new_variant["Last Edited"] <= old_variant["Last Edited"]:
                    # if date_of_upload <= old_variant["Last Edited"]: # uncomment this line to FORCE UPDATE of all variants
                    up_to_date = True
                    break
                else:
                    if old_variant["SCV"] == "": break #missing SCV for them, will likely result in error when uploaded
                    new_variant["SCV"] = old_variant["SCV"]
                    new_variant["Last Edited"] = date_of_upload
                    variants_update.append(new_variant)
                    to_update = True
                    break
        if not to_update and not up_to_date:  # it's a new variant (no SCV)
            variants_upload.append(new_variant)
    return variants_upload, variants_update
