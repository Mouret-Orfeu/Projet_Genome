from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import reverse_complement, Seq, complement
from Bio.SeqFeature import CompoundLocation
import re
def downloadNCRegion(espece: str):

    Entrez.email = "vvxStcd42tuYvQD2B@gmail.com"
    search_term = f"{espece}[Organism] AND NC_000000:NC_999999[Accession]"

    # Effectuer la recherche dans la base de données NCBI
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax= 10 ** 9)
    record = Entrez.read(handle)

    # Récupérer les identifiants des fichiers NC trouvés
    id_list = record["IdList"]

    nb_NC = len(id_list)
    print(nb_NC)


    handle.close()

#downloadNCRegion("Bufo bufo")

def download_gb_file(accession_id: str, file_path: str):
    """
    Télécharge un fichier NC de GenBank à partir d'un ID d'accès et l'enregistre localement
    :param accession_id: l'ID d'accès pour le fichier GenBank à télécharger
    :param file_path: le chemin de fichier local pour enregistrer le fichier téléchargé
    """
    Entrez.email = "votre_email@example.com"  # Entrez exige une adresse e-mail pour accéder aux données
    with Entrez.efetch(db="nucleotide", id=accession_id, rettype="gbwithparts", retmode="text") as handle:
        record = SeqIO.read(handle, "genbank")
        with open(file_path, "w", encoding='utf-8') as output_file:
            SeqIO.write(record, output_file, "genbank")
    print(f"Fichier GenBank enregistré à : {file_path}")


#download_gb_file("NC_001371.1", "Ecoli_K12.gb")

def extract_sequences_cds(record):#record,espece_path, region, espece
    sequences = {}
    for feature in record.features:
        if feature.type == "CDS":
            if feature.location_operator == "join":
                parts =feature.location.parts
                parts_str_list = []
                for part in parts:
                    part_str = f"{part.start + 1}..{part.end}"
                    parts_str_list.append(part_str)
                parts_str = "join(" + ",".join(parts_str_list) + ")"
                seq_parts = []
                exon_num=0
                for part in feature.location.parts:
                    start_str = str(part.start)
                    end_str = str(part.end)
                    if re.match("^[0-9]+$", start_str) and re.match("^[0-9]+$", end_str):
                        start=int(start_str)
                        #print(start)
                        end=int(end_str)
                        #print(end)
                        strand = part.strand
                        seq_part = feature.extract(record).seq
                        #print(len(seq_part))
                        if len(seq_part)!=0:
                            #print(len(seq_part))
                            if start < end and start >=1:
                                if strand == -1:
                                   seq_part = reverse_complement(seq_part)
                                else:
                                    seq_part=seq_part
                            bounds = f"{parts_str}:{start+1}..{end} Exon {exon_num +1}" if strand == 1 else f"complement({parts_str}):complement({start+1}..{end}) Exon {exon_num +1}"
                            sequences[bounds] = seq_part
                            exon_num=exon_num+1

                                #write the sequence in a file .txt

                        seq_parts.append(seq_part)
                        #print(seq_parts)
                if(len(seq_parts)!=0):
                    cds_seq = Seq("")
                    #print(seq_parts)
                    #print(len(seq_parts))
                    for seq_part in seq_parts:
                        cds_seq += seq_part
                        #print(len(cds_seq))
                    if len(cds_seq) % 3 == 0: # remove last comma
                        bounds = f"{parts_str}" if feature.location.strand == 1 else f"complement({parts_str})"
                        sequences[bounds] = cds_seq

                    #write the sequence in a file .txt


            else:

                start_str = str(feature.location.start)
                end_str = str(feature.location.end)
                if re.match("^[0-9]+$", start_str) and re.match("^[0-9]+$", end_str):
                    start = int(start_str)
                    #print(start)
                    end = int(end_str)
                    strand = feature.location.strand
                    cds_seq = feature.extract(record).seq
                    #print(len(cds_seq))
                    if len(cds_seq) % 3==0:
                       # print(len(cds_seq))
                        if strand == -1:
                            if start < end and start >= 1 and end <= len(cds_seq):
                                    cds_seq = reverse_complement(cds_seq)
                        else:
                            cds_seq=cds_seq
                        bounds = f"{start+1}..{end}" if feature.location.strand == 1 else f"complement({start+1}..{end})"

                        sequences[bounds] = cds_seq

    return sequences



def extract_sequences_intron(record):

    sequences = {}
    for feature in record.features:
        if feature.type == "CDS":
            location_operator = feature.location_operator
            # Si l'opérateur join
            if location_operator == "join":
                parts = sorted(feature.location.parts, key=lambda p: p.start)
                #print(parts)
                parts_str_list = []
                for i in range(len(parts)):
                    part_str = f"{parts[i].start + 1}..{parts[i].end}"
                    parts_str_list.append(part_str)
                parts_str = "join(" + ",".join(parts_str_list) + ")"
                for i in range(1, len(parts)):
                    # Position de fin de l'exon précédent
                    prev_end_str = str(parts[i-1].end)
                    cur_star_str =str(parts[i].start)
                    #print(prev_end)
                    if re.match("^[0-9]+$", prev_end_str) and re.match("^[0-9]+$", cur_star_str):
                        prev_end=int(prev_end_str)
                        cur_start=int(cur_star_str) +1
                        # Position de début de l'exon courant
                        #print(cur_start)
                       # Positions de l'intron
                        intron_start = prev_end + 1
                        intron_end = cur_start - 1
                        print(intron_start,intron_end)

                        intron_seq=feature.extract(record).seq[intron_start:intron_end]
                        if(intron_start< intron_end) and prev_end < cur_start:
                            if parts[i-1].strand == -1:
                                intron_seq =reverse_complement(intron_seq)
                            else:
                                intron_seq=intron_seq
                            print(len(intron_seq))

                            if len(intron_seq)%3==0:# Séquence de l'intron
                                if parts[i-1].strand == -1:
                                   bounds = f"complement({parts_str}) :Intron {i} complement({str(int(intron_start))}" + str("..") + f"{intron_end})"
                                else:
                                   bounds = f"{parts_str} :Intron {i} {intron_start}..{intron_end}"
                                sequences[bounds] = intron_seq


    return sequences


#download_gb_file("NC_000962.3", "variole.gb")
#download_gb_file("NC_001802.1", "mosaic_tobacco_virus.gb")


# Ouvrir un fichier GenBank contenant plusieurs enregistrements
filename = "Ecoli_K12.gb"
record = SeqIO.read(filename, "genbank")

# Extraire les séquences CDS
sequences_cds = extract_sequences_cds(record)
print(sequences_cds)

# Afficher les séquences CDS extraites
#for key, value in sequences_cds.items():
 #  print(key, value)##
def downloadNCRegion(espece):
    search_term = f" {espece}[Organism]) AND NC_000000:NC_999999[Accession]"

    Entrez.email = "vvxStcd42tuYvQD2B@gmail.com"



    # Effectuer la recherche dans la base de données NCBI
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10**9)
    record = Entrez.read(handle)

    # Récupérer les identifiants des fichiers NC trouvés
    id_list = record["IdList"]

    nb_NC = len(id_list)

    handle.close()

    # Pour chaque identifiant, récupérer les données et stocker en RAM avec StringIO
    records = []
    for i, record_id in enumerate(id_list):
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gbwithparts", retmode="text")
        record = SeqIO.read(StringIO(handle.read()), "genbank")
        handle.close()
        print(f"{record.id}\n")

            # Extraire les séquences de la région de l'enregistrement en appelant la fonction correspondante (ça les écrit dans un fichier aussi)
        sequences=extract_sequences_cds(record)

        records.append(sequences)

        # Écrire les séquences dans un fichier texte
        for seq_name, seq in sequences.items():
              header = f"{seq_name}"
              filename = f"CDS_{espece}_{record.id}.txt"
              with open(filename, "a", encoding='utf-8') as f:
                f.write(f"{header}\n{seq}\n")

    print("fin downloadNC")
    return
#downloadNCRegion("Escherichia coli")
