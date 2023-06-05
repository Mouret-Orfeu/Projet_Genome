import os
from typing import List

from gui import MainDialog
from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import  Seq,reverse_complement
import re

#extraire les CDS
def extract_sequences_cds(record, organism_path_list: List[str], region: str, organism: str, main_dialog: MainDialog):#record,espece_path, region, espece

    #if file " filename = f"Results/{espece_path}/{region}_{espece}_{record.id}.txt" already exists, return
    os.path.join
    filename = os.path.join("Results", *organism_path_list, f"{region}_{organism}_{record.id}.txt")
    
    #log filemname
    #main_dialog.log_signal.emit(f"filename: {filename}", {"color":"lightgreen"})

    if os.path.exists(filename):
        main_dialog.log_signal.emit(f"Skipped {region} sequences (already treated)", {"color":"lightgreen"})
        return

    main_dialog.log_signal.emit(f"Extracting {region} sequences...", {"color":"lightgreen"})

    #nombre de sequence de CDS trouvées
    nb_CDS_seq=0

    sequences = {}
    for feature in record.features:
        if feature.type == "CDS":
            if feature.location_operator == "join":
                parts = feature.location.parts
                parts_str_list = []
                for part in parts:
                    part_str = f"{part.start + 1}..{part.end}"
                    parts_str_list.append(part_str)
                parts_str = "join(" + ",".join(parts_str_list) + ")"
                seq_parts = []
                prev_end = 0
                exon_num=0
                for part in parts:
                    start_str = str(part.start)
                    end_str = str(part.end)
                    if re.match("^[0-9]+$", start_str) and re.match("^[0-9]+$", end_str):
                        start=int(start_str)
                        end=int(end_str)
                        strand = part.strand
                        seq_part = record.seq[start:end]
                        #print(len(seq_part))
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
                        seq_name= bounds
                        header = f"CDS {organism} {record.id}:{seq_name}"
                        write_seq_in_file(filename, header, seq_part)

                        seq_parts.append(seq_part)
                if(len(seq_parts)!=0):            
                    cds_seq = Seq("")
                    for seq_part in seq_parts:
                        cds_seq += seq_part
                    if len(cds_seq) % 3 == 0: # remove last comma
                        bounds = f"{parts_str}" if feature.location.strand == 1 else f"complement({parts_str})"
                        sequences[bounds] = cds_seq

                        #write the sequence in a file .txt
                        seq_name= bounds
                        header = f"CDS {organism} {record.id}:{seq_name}"
                        write_seq_in_file(filename, header, cds_seq)

            else:
                start_str = str(feature.location.start)
                end_str = str(feature.location.end)
                if re.match("^[0-9]+$", start_str) and re.match("^[0-9]+$", end_str):
                    start = int(start_str)
                    end = int(end_str)
                    strand = feature.location.strand
                    cds_seq = feature.extract(record).seq
                    if len(cds_seq) % 3==0:
                        if strand == -1:
                            if start < end and start >= 1 and end <= len(cds_seq):
                                cds_seq = reverse_complement(cds_seq)
                        else:
                            if start < end and start >= 1 and end <= len(cds_seq):
                                cds_seq=cds_seq    
                        bounds = f"{start+1}..{end}" if feature.location.strand == 1 else f"complement({start+1}..{end})"

                        sequences[bounds] = cds_seq
                        seq_name = bounds
                        header = f"CDS {organism} {record.id}:{seq_name}"
                        write_seq_in_file(filename, header, cds_seq)
                        nb_CDS_seq+=1
    if nb_CDS_seq==0:
        main_dialog.log_signal.emit("No CDS sequence found", {"color":"red"})


    return


#extraire les introns
def extract_sequences_intron(record, organism_path_list: List[str], region: str, organism: str, main_dialog: MainDialog):

    #if file " filename = f"Results/{espece_path}/{region}_{espece}_{record.id}.txt" already exists, return
    filename = os.path.join("Results", *organism_path_list, f"{region}_{organism}_{record.id}.txt")

    if os.path.exists(filename):
        main_dialog.log_signal.emit(f"Skipped {region} sequences (already treated)", {"color":"lightgreen"})
        return

    main_dialog.log_signal.emit(f"Extracting {region} sequences...", {"color":"lightgreen"})


    sequences = {}

    #conteur du nombre d'intron
    nb_intron_seq=0

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
                intron_num=0
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

                        intron_seq=record.seq[intron_start:intron_end]
                        if(intron_start< intron_end) and prev_end < cur_start:
                            if parts[i-1].strand == -1:
                                intron_seq =reverse_complement(intron_seq)
                            else:
                                intron_seq=intron_seq
                            print(len(intron_seq))

                            # Séquence de l'intron
                            if parts[i-1].strand == -1:
                                bounds = f"complement({parts_str}) :Intron {intron_num+1} complement({str(int(intron_start))}" + str("..") + f"{intron_end})"
                            else:
                                bounds = f"{parts_str} :Intron {intron_num+1} {intron_start}..{intron_end}"
                            sequences[bounds] = intron_seq
                            intron_num=intron_num + 1

                            #write the sequence in a file .txt
                            seq_name= bounds
                            header = f"intron {organism} {record.id}:{seq_name}"
                            write_seq_in_file(filename, header, intron_seq)
                            nb_intron_seq+=1
    if nb_intron_seq==0:
        main_dialog.log_signal.emit("No intron found", {"color":"lightgreen"})
    return


#download Fichiers NC Ram bon
def downloadNCRam(organism: str):
    # Spécifier votre adresse email pour vous connecter à GenBank
    Entrez.email = "votre_email@votre_domaine.com"

    # Créer la chaîne de recherche
    search_term = organism + "[Organism] AND (NC_*[Accession])"

    # Effectuer la recherche dans la base de données NCBI
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10 ** 9)
    record = Entrez.read(handle)

    # Récupérer les identifiants des fichiers NC trouvés
    id_list = record["IdList"]
    handle.close()

    # Pour chaque identifiant, récupérer les données et stocker en RAM avec StringIO
    records = []
    for i, record_id in enumerate(id_list):
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gbwithparts", retmode="text")
        record = SeqIO.read(StringIO(handle.read()), "genbank")
        records.append(record)
        handle.close()
        print(f"{record.id}\n")
    print(f"{len(records)} fichiers NC ont été téléchargés et stockés en RAM.")
    return records

#records = downloadNCRam("Escherichia coli")#,espece_path,region,espece


def download(organism: str):
    # Spécifier votre adresse email pour vous connecter à GenBank
    Entrez.email = "mariam.diakite@etu.unistra.fr"

    # Créer la chaîne de recherche
    search_term = organism + "[Organism] AND (NC_*[Accession])"

    # Effectuer la recherche dans la base de données NCBI
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10 ** 9)
    record = Entrez.read(handle)

    # Récupérer les identifiants des fichiers NC trouvés
    id_list = record["IdList"]
    handle.close()

    # Pour chaque identifiant, récupérer les données et stocker en RAM avec StringIO
    records = []
    for i, record_id in enumerate(id_list):
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gbwithparts", retmode="text")
        record = SeqIO.read(StringIO(handle.read()), "genbank")
        handle.close()
        # Extraire les séquences CDS complètes de l'enregistrement et stocker dans un dictionnaire
        sequences = extract_sequences_cds(record)
        print(sequences)
        # Ajouter le dictionnaire des séquences CDS à la liste des enregistrements
        records.append(sequences)

        # Écrire les séquences CDS dans un fichier texte
        for cds_name, cds_seq in sequences.items():
            cds_header = f"CDS {organism} {record.id}: {cds_name}"
            filename = f"CDS_{organism}_{record.id}.txt"
            with open(filename, "a", encoding='utf-8') as f:
                f.write(f"{cds_header}\n{cds_seq}\n")
        print( f"{record.id}\n")
    print(f"{len(records)} fichiers NC ont été téléchargés et les séquences CDS ont été écrites dans les fichiers textes.")

    return records

#records = download("Escherichia coli")

extract_functions = {
    "CDS": extract_sequences_cds,
    "intron": extract_sequences_intron,
    # Ajouter d'autres fonctions d'extraction pour d'autres régions si nécessaire
}
def downloadNCRegion(organism_path_list: List[str], regions: str, main_dialog: MainDialog):

    organism = organism_path_list[-1]
    main_dialog.log_signal.emit(f"Organism: <em style=\"color: lightblue;\">{organism}</em>", {"bold":True})

    progress_value = 0

    temp = " AND ".join(f"({s}[Organism])" for s in organism_path_list)
    search_term = f"{temp} AND NC_000000:NC_999999[Accession]"

    Entrez.email = "vvxStcd42tuYvQD2B@gmail.com"

    if main_dialog.processing_thread.stop_requested:
        return

    # Effectuer la recherche dans la base de données NCBI
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10 ** 9)
    record = Entrez.read(handle)

    # Récupérer les identifiants des fichiers NC trouvés
    id_list = record["IdList"]

    nb_NC = len(id_list)
    nb_regions = len(regions)
    steps = nb_NC*nb_regions
    main_dialog.set_max_progress_bar2_signal.emit(100)


    #if no NC found return and print in main dialog
    if not id_list:
        main_dialog.log_signal.emit("No NC found", {"color":"lightgreen"})
        return

    handle.close()

    main_dialog.set_value_progress_bar2_signal.emit(0)

    # Pour chaque identifiant, récupérer les données et stocker en RAM avec StringIO
    records = []
    n = len(id_list)
    for i, record_id in enumerate(id_list, start=1):
        if main_dialog.processing_thread.stop_requested:
            return

        main_dialog.log_signal.emit(f"NC file {i}/{n}", {"bold":True})

        main_dialog.log_signal.emit("Downloading NC file", {"color":"lightgreen"})

        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gbwithparts", retmode="text")
        record = SeqIO.read(StringIO(handle.read()), "genbank")
        print(f"{record.id}\n")

        # #print NC in main dialog
        # main_dialog.log_signal.emit(f"Extracting sequences from {record.id} (file {i}/{n})...", {"color":"lightgreen"})

        for region in regions:
            if main_dialog.processing_thread.stop_requested:
                return

            if region in extract_functions.keys():
                # Extract sequences by calling the appropriate function and write them in a file
                extract_functions[region](record, organism_path_list, region, organism, main_dialog)
            else:
                main_dialog.log_signal.emit(f"Skipping {region} sequences (not implemented)", {"color":"orange"})

            progress_value += 1
            main_dialog.set_value_progress_bar2_signal.emit(progress_value*100//steps)

    handle.close()
    print("fin downloadNCRegion")
    #print NC in main dialog


def write_seq_in_file(filename: str, header: str, seq: str):
    with open(filename, "a", encoding='utf-8') as f:
        f.write(f">{header}\n{seq}\n")


#if __name__ == "__main__":
#    records = download("Escherichia coli")
