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

def compare_variants(new_variants, old_variants, date_of_upload):
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
            old_variant_coords = "\t".join(
                [
                    old_variant["#Chromosome"],
                    old_variant["Start"],
                    old_variant["Stop"],
                    old_variant["Ref/Alt"],
                ]
            )
            if (new_variant["hgvs c."] == old_variant["hgvs c."] != "") or (
                new_variant_coords == old_variant_coords
            ):
                if new_variant["Last Edited"] <= old_variant["Last Edited"]:
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

def compare_haplotypes(new_haplos, old_haplos):
    """Function to check whether variants are to be uploaded as novel haplotypes or updated in ClinVar.
    The function checks if 
    
    Parameters:
        new_haplos: list containing haplotypes from the input file.
        old_haplos: list containing haplotypes from the reference file.
    
    Output:
        haplos_novel: list containing haplotypes to be uploaded as novel.
        haplos_update: list containing haplotypes to be updated."""
    
    old_haplos = read_file(old_haplos)
    haplos_novel = {}
    haplos_update = {}
    for new_haplo in new_haplos:
        to_update = False
        up_to_date = False
        new_associated_variants, new_haplo_hgvsc, new_haplo_hgvsp, new_haplo_classification, new_upload_type, new_last_edited = new_haplos[new_haplo].strip(" ").split("; ")
        for old_haplo in old_haplos:
            if new_haplo == old_haplo["hgvs c."]: #new_haplo is hgvs c., key of the haplotype dictionary
                old_haplo_hgvsc, old_haplo_classification, old_associated_variants, old_haplo_hgvsp, old_last_edited, old_haplo_SCV = old_haplo.values()
                if new_last_edited <= old_last_edited: #up to date
                    up_to_date = True
                    break
                else: 
                    new_haplo_info = "; ".join([new_haplos[new_haplo], old_haplo_SCV])
                    to_update = True
                    haplos_update[new_haplo] = new_haplo_info
                    break
        if not to_update and not up_to_date: # new haplotype
            haplos_novel[new_haplo] = new_haplos[new_haplo]
    return haplos_novel, haplos_update
