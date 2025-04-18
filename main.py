from Bio import Entrez
from Bio import SeqIO
import xmltodict
import time
import configparser

conf = configparser.ConfigParser()
conf.read("conf.ini")
Entrez.email = conf.get("General", "email")

# Variables
snp_dict = {}
snp_skipped = []
gene_list = []
cs_list = ["PATHOGENIC", "ASSOCIATION", "OTHER", "UNCERTAIN", "BENIGN", "NONE"]
clinical_signifiance_dict = {'pathogenic':cs_list[0], 'pathogenic-likely-pathogenic':cs_list[0], 'likely-pathogenic':cs_list[0], 
'pathogenic-low-penetrance':cs_list[0], 'affects':cs_list[1], 'risk-factor':cs_list[1], 'drug-response':cs_list[1], 
'association':cs_list[1], 'other':cs_list[2], 'uncertain-significance':cs_list[3], 'conflicting-interpretations-of-pathogenicity':cs_list[3], 
'no-classifications-from-unflagged-records':cs_list[3], 'not-provided':cs_list[3], 'likely-benign':cs_list[3], 
'benign-likely-benign':cs_list[3], 'benign':cs_list[3], 'None':cs_list[4]}


with open("input.txt", "r", encoding="utf-8") as file:
    lines = file.read().splitlines()
    for line in lines:
        snp_id = line.split()[0]
        match snp_id[0]:
            case "#":
                continue
            case "i":
                genotype = line.split()[3]
                snp_skipped.append({"snp":snp_id, "genotype":genotype})
                continue
            case _:
                genotype = line.split()[3]
                snp_dict[snp_id[2:]] = genotype


def fetch_snp(id_dict):
    search_results = Entrez.read(Entrez.epost("snp", id=",".join(id_dict.keys())))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]
    
    batch_size = 100000
    count = len(id_dict.keys())
    result = []
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Downloading snp %i to %i" % (start + 1, end))
        stream = Entrez.efetch(
            db="snp", 
            webenv=webenv, 
            query_key=query_key, 
            rettype="gb", 
            retmode="text", 
            retstart=start,
            retmax=batch_size,
            idtype="acc",
        )
        data = stream.read()
        stream.close()
        print("Processing SNPs")
        result.extend(extract_snp(data)) #Parse using Entrez
    return result

def extract_snp(data):
    top = xmltodict.parse(f"<top>{data}</top>")['top']
    res = []
    # Locate SNPs
    for snp in top["DocumentSummary"]:
        entry = {}
        entry["SNP_ID"] = snp.get("SNP_ID")
        if entry["SNP_ID"] == None:
            continue
        entry["CLINICAL_SIGNIFICANCE"] = get_clinical_signifiance_category(snp["CLINICAL_SIGNIFICANCE"])
        entry["SPDI"] = snp["SPDI"]
        entry["GENES"] = snp["GENES"]
        # Get a record of every gene for requesting detail later
        genes = entry["GENES"]
        if genes is not None:
            genes = genes["GENE_E"]
            if not isinstance(genes, list):
                genes = [genes]
            for gene in genes:
                gene_list.append(gene["GENE_ID"])
        res.append(entry)

    return res

def fetch_genes(gene_set):
    print(f"Fetching {len(gene_set)} genes, be patient the process can be very long")
    search_results = Entrez.read(Entrez.epost("gene", id=",".join(gene_set)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]
    
    handle = Entrez.efetch(
            db="gene",
            id=gene_set,
            retmode="xml"
        )
    return Entrez.read(handle, validate=False)

# creates an array of all the genes with snp info
def add_genes(snp_data, gene_data):
    print(f"Fusing {len(gene_data)} genes with SNPs")
    gene_array = []
    gene_array_mutated = []
    steps = 100
    for i in range(len(gene_data)):
        gene = gene_data[i]
        if i % steps == 0:
            percentage = "{:.2f}".format(i / len(gene_data) * 100)
            print(f"Fusion progress {percentage} %")
        if "Entrezgene_prot" in gene:
            protein_desc = gene["Entrezgene_prot"]["Prot-ref"].get("Prot-ref_desc") or\
                ",".join(gene["Entrezgene_prot"]["Prot-ref"].get("Prot-ref_name"))
        summary = gene.get("Entrezgene_summary") or ""
        gene_desc = gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_desc"]
        gene_locus = gene["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
        gene_id = gene["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]

        snp_clinical_signifiance = {k: [] for k in cs_list}

        snp_clinical_signifiance_mutated = {k: [] for k in cs_list}
        mutated = False
        for snp in snp_data:
            genes_in_snp = snp["GENES"]
            if genes_in_snp is None:
                continue
            genes_content = genes_in_snp["GENE_E"]
            if not isinstance(genes_content, list):
                genes_content = [genes_content]
            for gene_content in genes_content:
                snp_gene_id = gene_content["GENE_ID"]
                genotype = snp_dict.get(snp["SNP_ID"])
                if gene_id == snp_gene_id and genotype is not None:
                    reference_genotype = snp["SPDI"].split(":")[2]
                    snp_mutated = False
                    # Check genotype against reference 
                    for char in genotype:
                        snp_mutated = snp_mutated or char != reference_genotype
                    if snp_mutated:
                        mutated = True
                        snp_clinical_signifiance_mutated[snp["CLINICAL_SIGNIFICANCE"]]\
                        .append(f"{snp["SNP_ID"]} {genotype}")
                    snp_clinical_signifiance[snp["CLINICAL_SIGNIFICANCE"]]\
                    .append(f"{snp["SNP_ID"]} {genotype}:{"Mutated" if mutated else "Reference"}")
                    break

        gene_base = {"gene_id":gene_id, "gene_locus":gene_locus, "gene_desc":gene_desc, "summary":summary.replace('"','""')}
        gene_info = gene_base.copy()
        for col in cs_list:
            gene_info[col] = snp_clinical_signifiance[col]
        gene_array.append(gene_info)
        if(mutated):
            gene_info_mutated = gene_base.copy()
            for col in cs_list:
                gene_info_mutated[col] = snp_clinical_signifiance_mutated[col]
            gene_array_mutated.append(gene_info_mutated)
    return [gene_array, gene_array_mutated]

# Find the most relevant clinical signifiance and returns a category for it
def get_clinical_signifiance_category(cs_string):
    cs_string = cs_string or "None"
    cs_list = cs_string.split(",")
    best = ""
    best_pos = 999
    for cs in cs_list:
        pos = 0
        for key in clinical_signifiance_dict.keys():
            pos = pos + 1
            if key == cs and pos < best_pos:
                best = clinical_signifiance_dict[key]
                best_pos = pos
                break
    return best

def export_csv(data, filename = "output.csv"):
    # Header
    csv_str = ",".join(str(x) for x in data[0].keys()) + "\n"
    # Content
    for element in data:
        csv_str = f"{csv_str}{",".join(f"\"{x}\"" for x in element.values())}\n"
    f = open(filename, "a")
    f.write(csv_str)
    f.close()

start = time.time()
snp_data = fetch_snp(snp_dict)

final_data = add_genes(snp_data, fetch_genes(set(gene_list)))
print("Exporting output.csv")
export_csv(final_data[0])
export_csv(final_data[1], "output_mutated.csv")
export_csv(snp_skipped, "skipped.csv")

end = time.time()
print(f"Done in {end - start} seconds")
