import urllib.request
import os
import csv
import shutil
from gui import MainDialog

def read_tsv_file(file_path):
    data = []
    with open(file_path, "r", encoding='utf-8') as f:
        reader = csv.reader(f, delimiter="\t")
        next(reader)  # skip header row
        for row in reader:
            try:
                organism, kingdom, group, subgroup, *_ = row
            except ValueError:
                # Skip rows that don't contain the expected number of columns
                continue
            organism = organism.replace(':', '').replace('/', '').replace('?', '')
            data.append((kingdom, group, subgroup, organism))
    return data

def create_arborescence(main_dialog: MainDialog):
    # Path to the TSV file
    overview = "data/overview.txt"
    new_overview = "data/overview_new.txt"

    # download the TSV file and save it in data/overview_new.txt
    main_dialog.log_signal.emit("Downloading `overview.txt`...", {"color":"lightgreen"})
    os.makedirs("data", exist_ok=True)
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt'
    urllib.request.urlretrieve(url, new_overview)

    if os.path.exists(overview):
        main_dialog.log_signal.emit("Updating arborescence...", {"color":"lightgreen"})

        old_data = set(read_tsv_file(overview))
        new_data = set(read_tsv_file(new_overview))

        added_data = new_data - old_data
        removed_data = old_data - new_data

        for kingdom, group, subgroup, organism in removed_data:
            path = os.path.join("Results", kingdom, group, subgroup, organism)
            if os.path.exists(path):
                shutil.rmtree(path)

        for kingdom, group, subgroup, organism in added_data:
            path = os.path.join("Results", kingdom, group, subgroup, organism)
            os.makedirs(path, exist_ok=True)

        os.replace(new_overview, overview)
    else:
        main_dialog.log_signal.emit("Creating arborescence...", {"color":"lightgreen"})

        data = read_tsv_file(new_overview)

        for kingdom, group, subgroup, organism in data:
            path = os.path.join("Results", kingdom, group, subgroup, organism)
            os.makedirs(path, exist_ok=True)

        os.replace(new_overview, overview)
