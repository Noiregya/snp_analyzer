from Bio import Entrez
from Bio import SeqIO
import xmltodict
import time
import configparser

conf = configparser.ConfigParser()
conf.read("conf.ini")
Entrez.email = conf.get("General", "email")

snp_ids = []
snp_id_skipped = []

with open("input.txt", "r", encoding="utf-8") as file:
    lines = file.read().splitlines()
    for line in lines:
        snp_id = line.strip().split()[0]
        match snp_id[0]:
            case "#":
                continue
            case "i":
                snp_id_skipped.append(snp_id)
                continue
            case _:
                snp_ids.append(snp_id[2:])

dataset = []

def fetch_snp(id_list):
    search_results = Entrez.read(Entrez.epost("snp", id=",".join(id_list)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]
    
    batch_size = 100000
    count = len(id_list)
    for start in range(0, count, batch_size):
        end = min(count, start + batch_size)
        print("Going to download record %i to %i" % (start + 1, end))
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
        extract(data)

#def fetch_gene():

def extract(data):
    top = xmltodict.parse(f"<top>{data}</top>")['top']
    for snp in top["DocumentSummary"]:
        break #TO REMOVE
        genes = snp["GENES"]
        if genes is not None:
            g_names = ''
            for gene in genes.keys():
                g_names = genes[gene]
            #print(snp["GENES"])

start = time.time()
#fetch_snp(snp_ids)
fetch_snp(["548049170",
"12210761",
"11755766",
"1029124",
"4712416",
"6926746",
"34825161",
"17631613",
"75072109",
"12201803",
"77217387",
"9368113",
"6899608",
"56949706",
"74400123",
"115798487",
"10806909",
"1413188",
"115267294",
"7452838",
"75514178",
"12111456",
"6900553",
"4712449",
"4712456",
"6456301",
"6456302",
"16893173",
"12191126",
"9350235",
"11963800",
"6940110",
"56218989",
"79541139",
"1333652",
"9465765",
"6456357",
"12525060"])
end = time.time()
print(f"done in {end - start} seconds")

#def get_accessions(genetic_accession):
#    handle = Entrez.efetch(db="snp", id=genetic_accession,rettype="gb", retmode="text")
#    record = SeqIO.read(handle, "genbank")
#    dbsource = record.annotations.get("db_source")
#    dbsource = dbsource.replace("accession ", "")
#    dbsource = dbsource.replace("embl ", "")
#    locus_tag = record.features[-1].qualifiers.get("locus_tag", [""])[0]
#    if locus_tag == "":
#         locus_tag = record.features[-1].qualifiers.get("gene", [""])[0]
#     handle.close()
#     return dbsource, locus_tag

# for genetic_accession in genetic_accessions:
#     nuc_accession, locus_tag = get_accessions(genetic_accession)
#     ID_delimited = f"{genetic_accession} {nuc_accession} {locus_tag}"
#     accession_all.append(ID_delimited)

# output_file = open("D:/Bioinformatics/nucleotide_accession.txt", "w", encoding="utf-8")
# for x in accession_all:
#     output_file.write(x + "\n")
# output_file.close()

