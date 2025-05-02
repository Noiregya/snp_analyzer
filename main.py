import http
import time
import json
import os
import io
import configparser
import xmltodict
from Bio import Entrez
from Bio import SeqIO

conf = configparser.ConfigParser()
conf.read("conf.ini")
Entrez.email = conf.get("General", "email")
feature_types = ["gene"]

# Variables
snp_skipped = []
gene_list = []
cs_list = ["PATHOGENIC", "ASSOCIATION", "OTHER", "UNCERTAIN", "BENIGN", "NONE"]
clinical_signifiance_dict = {'pathogenic':0, 'pathogenic-likely-pathogenic':0, 'likely-pathogenic':0, 
'pathogenic-low-penetrance':0, 'affects':1, 'risk-factor':1, 'drug-response':1, 
'association':1, 'other':2, 'uncertain-significance':2, 'conflicting-interpretations-of-pathogenicity':2, 
'no-classifications-from-unflagged-records':2, 'not-provided':2, 'likely-benign':3, 
'benign-likely-benign':3, 'benign':3, 'None':4}

assembly_grch37_chromosomes_dict = {"1":["CM000663.1","NC_000001.10"],
"2":["CM000664.1","NC_000002.11"],
"3":["CM000665.1","NC_000003.11"],
"4":["CM000666.1","NC_000004.11"],
"5":["CM000667.1","NC_000005.9"],
"6":["CM000668.1","NC_000006.11"],
"7":["CM000669.1","NC_000007.13"],
"8":["CM000670.1","NC_000008.10"],
"9":["CM000671.1","NC_000009.11"],
"10":["CM000672.1","NC_000010.10"],
"11":["CM000673.1","NC_000011.9"],
"12":["CM000674.1","NC_000012.11"],
"13":["CM000675.1","NC_000013.10"],
"14":["CM000676.1","NC_000014.8"],
"15":["CM000677.1","NC_000015.9"],
"16":["CM000678.1","NC_000016.9"],
"17":["CM000679.1","NC_000017.10"],
"18":["CM000680.1","NC_000018.9"],
"19":["CM000681.1","NC_000019.9"],
"20":["CM000682.1","NC_000020.10"],
"21":["CM000683.1","NC_000021.8"],
"22":["CM000684.1","NC_000022.10"],
"X":["CM000685.1","NC_000023.10"],
"Y":["CM000686.1","NC_000024.9"],
"MT":["", "NC_012920.1"]}


def read_input(input_file):
    snp_elements = []
    with open(input_file, "r", encoding="utf-8") as file:
        lines = file.read().splitlines()
        for line in lines:
            splitline = line.split()
            if splitline[0][0] == "#" or splitline[0] == "rsid":
                continue
            rsid, chromosome, position, genotype = line.split()
            elem = {"rsid":rsid,"chromosome":chromosome,"position":position,"genotype":genotype}
            if genotype == "--":
                snp_skipped.append(elem)
                continue
            snp_elements.append(elem)
    return snp_elements


def get_chromosomes(chromosome_dict):
    # Download chromosomes the first time
    chromosome_dict = chromosome_dict.copy()
    if not os.path.isfile("chromosomes.json"):
        chromosomes = {}
        for key in chromosome_dict.keys():
            features = []
            print(f"Caching chromosome features... Chromosome {key}")
            with Entrez.efetch(db="nuccore", rettype="gbwithparts", retmode="text", id=chromosome_dict[key][1]) as handle:
                chromosome = SeqIO.read(handle, "gb")
                for feature in chromosome.features:
                    if feature.type in feature_types:
                        features.append({"location_start": int(feature.location.start), "location_end": int(feature.location.end), "qualifiers": feature.qualifiers, "type": feature.type})
            chromosomes[key] = {"features": features, "seq": str(chromosome.seq)}
        with io.open("chromosomes.json", "w", encoding="utf-8") as file:
            file.write(json.dumps(chromosomes))

    # Open chromosome file
    with open("chromosomes.json", "r", encoding="utf-8") as file:
        print("Loading chromosome features...")
        data = file.read()
        chromosomes = json.loads(data)
        for chromosome in chromosomes.keys():
            chromosome_dict[chromosome].append(chromosomes[chromosome])
    return chromosome_dict


def get_shallow_snp(snp_elements, chromosome_dict):
    shallow_snps = {}
    
    last_chromosome = ''
    last_index = 0
    
    for elem in snp_elements:
        rsid, chromosome, position, genotype = elem.values()
        position = int(position)
        chromosome_info = chromosome_dict[chromosome][2]["features"]
        seq = chromosome_dict[chromosome][2]["seq"]
        elem = {"rsid":rsid,"chromosome":chromosome,"position":position,"genotype":genotype,"reference":seq[int(position - 1)]}

        if chromosome != last_chromosome:
            last_index = 0
        i = last_index

        for i in range(last_index, len(chromosome_info)):
            feature = chromosome_info[i]
            if feature["type"] != "gene": #Feature isn't gene 
                continue #Check next feature
            if position < feature["location_start"]: #This nucleotide is not contained within a gene
                break #Skip to next nucleotide
            if position <= feature["location_end"]: #This nucleotide is contained in the current gene
                ql = feature["qualifiers"]
                gene_id = list(filter(lambda x: x.split(":")[0] == "GeneID", ql["db_xref"]))[0].split(":")[1]
                gene_locus = ql["gene"][0]
                note = ql.get("note")
                gene_desc = ""
                if note is not None:
                    gene_desc = note[0]
                indel = False
                for char in genotype:
                    if char == "D" or char == "I":
                        indel = True
                        break

                elem["rsid"]=rsid
                elem["chromosome"] = chromosome
                elem["position"] = position
                elem["genotype"] = genotype
                elem["gene_id"] = gene_id
                elem["gene_locus"] = gene_locus
                elem["gene_desc"] = gene_desc
                elem["indel"] = indel
               
                gene_list.append(gene_id)
                shallow_snps[elem["rsid"]] = elem
                break
        last_index = i
        last_chromosome = chromosome
    return shallow_snps


def fetch_snp(snp_dict):
    snp_ids = []
    for snp in snp_dict.keys():
        if snp[:2] == "rs":
            snp_ids.append(snp[2:])
    if not os.path.isfile("snps.xml"):
        search_results = Entrez.read(Entrez.epost("snp", id=",".join(snp_ids)))
        query_key = search_results["QueryKey"]
        webenv = search_results["WebEnv"]

        data = ""
        batch_size = 100000
        count = len(snp_ids)
        for snp_start in range(0, count, batch_size):
            snp_end = min(count, snp_start + batch_size)
            print(f"Caching snp {snp_start + 1} to {snp_end}")
            stream = Entrez.efetch(
                db="snp",
                webenv=webenv,
                query_key=query_key,
                rettype="gb",
                retmode="text",
                retstart=snp_start,
                retmax=batch_size,
                idtype="acc",
            )
            data = data + stream.read()
            stream.close()
        print("Saving SNPs")
        with io.open("snps.xml", "w", encoding="utf-8") as file:
            file.write(f"<top>{data}</top>")

    # Open snp file
    with open("snps.xml", "r", encoding="utf-8") as file:
        print("Loading snps...")
        data = file.read()

    return assemble_snp(data, snp_dict) #TODO: Parse using Entrez


def assemble_snp(data, snp_dict):
    top = xmltodict.parse(data)['top']
    del data
    # Locate SNPs
    for snp in top["DocumentSummary"]:
        snp_id = snp.get("SNP_ID")
        if snp_id is None:
            continue
        rsid = f"rs{snp_id}"
        dict_entry = snp_dict.get(rsid)
        if dict_entry is None:
            continue
        clinical_signifiance = get_clinical_signifiance_category(snp["CLINICAL_SIGNIFICANCE"])
        spdi = snp["SPDI"]
        reference_genotype = ""
        if spdi is not None:
            reference_genotype = snp["SPDI"].split(":")[2]

        if reference_genotype != snp_dict[rsid]["reference"] and snp_dict[rsid]["indel"] is False:
            print(f"Alignement error: snp: {rsid}, snp database:{reference_genotype}, from sequence:{snp_dict[rsid]['reference']}")
            snp_dict[rsid]["reference"] = reference_genotype

        snp_dict[rsid]["clinical_signifiance"] = clinical_signifiance
    del top
    return snp_dict


def fetch_genes(snp_dict, gene_set):
    gene_set = list(gene_set)
    print(f"Fetching {len(gene_set)} genes, be patient the process can be very long")
    search_results = Entrez.read(Entrez.epost("gene", id=",".join(gene_set)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]

    batch_size = 500
    count = len(gene_set)
    for gene_start in range(0, count, batch_size):
        gene_end = min(count, gene_start + batch_size)
        attempts = 5
        i = 0
        while i < attempts:
            try:
                print(f"{(gene_start / count * 100):.2f}% Downloading gene {gene_start + 1} to {gene_end}")
                stream = Entrez.efetch(
                    db="gene",
                    retmode="xml",
                    complexity="4",
                    webenv=webenv,
                    query_key=query_key,
                    retstart=gene_start,
                    retmax=batch_size,
                )
                data, mutated_data = merge_snp_genes(snp_dict, Entrez.read(stream))
            except (Entrez.Parser.ValidationError, http.client.IncompleteRead) as err:
                i = i + 1
                print(f"{err}\nDownload error, trying again...")
                continue
            except RuntimeError as err:
                print(f"{err}\nDownload error, trying again...")
                search_results = Entrez.read(Entrez.epost("gene", id=",".join(gene_set)))
                query_key = search_results["QueryKey"]
                webenv = search_results["WebEnv"]
                continue
            break
        stream.close()
        export_csv(data)
        export_csv(mutated_data, "output_mutated.csv")
        print("Processing genes")
    return


def merge_snp_genes(snp_dict, gene_data):
    print(f"Fusing genes with SNPs")
    gene_dict = {}
    gene_dict_mutated = {}
    steps = 100
    i = 0
    length = len(gene_data)
    for gene in gene_data:
        gene_desc = gene["Entrezgene_gene"]["Gene-ref"].get("Gene-ref_desc")
        gene_locus = gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
        gene_id = gene["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
        summary = gene.get("Entrezgene_summary") or ""

        # Add snp information
        snp_clinical_signifiance = {k: [] for k in cs_list}
        snp_clinical_signifiance_mutated = {k: [] for k in cs_list}
        for snp in snp_dict.values():
            if snp.get("gene_id") is None or gene_id != snp.get("gene_id"):
                continue
            if gene_dict.get(snp["gene_id"]) is None:
                if gene_desc is None:
                    gene_desc = snp["gene_desc"]
                gene_dict[snp["gene_id"]] = {"gene_id":gene_id, "gene_locus":gene_locus, "gene_desc":gene_desc,
                "summary": summary, "mutated": False}

            snp_mutated = snp["indel"]
            genotype = snp["genotype"]
            rsid = snp["rsid"]
            reference_genotype = snp["reference"]
            for char in genotype:
                snp_mutated = snp_mutated or char != reference_genotype
            
            clinical_signifiance = snp.get("clinical_signifiance")
            if clinical_signifiance is not None:
                if snp_mutated:
                    gene_dict[snp["gene_id"]]["mutated"] = True
                    snp_clinical_signifiance_mutated[clinical_signifiance]\
                        .append(f"{rsid} {genotype}")
                snp_clinical_signifiance[clinical_signifiance]\
                    .append(f"{rsid} {genotype}:{'Mutated' if snp_mutated else 'Reference'}")

        # Create the entry if it's not in any snp
        gene_elem = gene_dict.get(gene_id)
        if gene_elem is None:
            gene_dict[gene_id] = {"gene_id":gene_id, "gene_locus":gene_locus, "gene_desc":gene_desc,
            "summary": summary, "mutated": False}
        
        # Add clinical signifiance columns
        for col in cs_list: #Convert list into string and add to gene_elem
            gene_elem[col] = ','.join(snp_clinical_signifiance[col])
        
        gene_dict[gene_id] = gene_elem
        if(gene_elem["mutated"]):
            gene_elem_mutated = gene_elem.copy()
            del gene_elem_mutated["mutated"]
            for col in cs_list:
                gene_elem_mutated[col] = ','.join(snp_clinical_signifiance_mutated[col])
            gene_dict_mutated[gene_id] = gene_elem_mutated
        i = i+1
    return [gene_dict, gene_dict_mutated]

# Find the most relevant clinical signifiance and returns a category for it
def get_clinical_signifiance_category(cs_string):
    cs_string = cs_string or "None"
    current_cs_list = cs_string.split(",")
    best = "NONE"
    best_pos = 999
    for cs in current_cs_list:
        for key in clinical_signifiance_dict.keys():
            pos = clinical_signifiance_dict[key]
            if key == cs and pos < best_pos:
                best_pos = pos
                best = cs_list[pos]
    return best


def export_csv(data, filename = "output.csv"):
    array = []
    if type(data) == dict:
        data = list(data.values())
    csv_str = ""
    # Header
    if not os.path.isfile(filename):
        csv_str = ",".join(f"\"{x}\"" for x in data[0].keys()) + "\n"
    # Content
    for element in data:
        csv_str = csv_str + ','.join('"' + str(x).replace('"','""') + '"' for x in element.values()) + "\n"

    f = open(filename, "a")
    f.write(csv_str)
    f.close()


start = time.time()
snp_elements = read_input("input.txt")#TODO: remove
chromosome_dict = get_chromosomes(assembly_grch37_chromosomes_dict)
shallow_snp = get_shallow_snp(snp_elements, chromosome_dict)
del chromosome_dict
snp_data = fetch_snp(shallow_snp)
export_csv(snp_skipped, "skipped.csv")
gene_set = set(gene_list)
fetch_genes(snp_data, gene_set)

print("Exporting output.csv")

end = time.time()
print(f"Done in {end - start} seconds")
