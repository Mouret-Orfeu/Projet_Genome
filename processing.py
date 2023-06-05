import os.path

from gui import MainDialog
from create_arborescence import create_arborescence
from extract_sequences import downloadNCRegion
from extract_sequences import extract_sequences_intron
from extract_sequences import extract_sequences_cds
from extract_sequences import extract_functions
from PyQt6.QtCore import QEventLoop


def setup(main_dialog: MainDialog):
    """This function is executed before the user can interact with the UI."""
    # Create an example arborescence. In the future,
    # the arborescence will be fetched from the website.
    main_dialog.log_signal.emit("Initial setup...", {"color":"green", "bold":True})

    # create_example_arborescence(main_dialog)
    create_arborescence(main_dialog)

    # Emit the signal to update the tree view
    main_dialog.updateTreeView_signal.emit()

    # Wait for the updateTreeView function to finish
    loop = QEventLoop()
    main_dialog.updateTreeView_done_signal.connect(loop.quit)
    loop.exec_()

    main_dialog.log_signal.emit("Ready and waiting.", {"color":"green", "bold":True})
    main_dialog.log_raw_signal.emit("")

def on_start_button_pressed(main_dialog: MainDialog):
    """This function is executed when the button Start is pressed."""
    main_dialog.log_signal.emit("Starting...", {"color":"green", "bold":True})

    main_dialog.set_value_progress_bar1_signal.emit(0)
    main_dialog.set_value_progress_bar2_signal.emit(0)
    main_dialog.set_format_progress_bar1_signal.emit("%v out of ??")

    #Get the selected regions
    selected_regions = [checkBox.text() for checkBox in main_dialog.functional_regions if checkBox.isChecked()]
    supplementary_region = main_dialog.textInput.text()
    if supplementary_region != "":
        selected_regions.append(supplementary_region)
    # main_dialog.log_signal.emit(f"[DEBUG] List of checked checkboxes: {selected_regions}", {"color": "orange"})

    # main_dialog.log_signal.emit(f'[DEBUG] Text input: "{main_dialog.textInput.text()}"', {"color": "orange"})

    #*********************************************************************************
    # Get the path of the selected file in the treeView widget
    # (suppose it is a file under the "Results" directory, and there is only 1 "Results" directory in the tree)

    #use the currentIndex method of the treeView to get the index of the currently selected item
    index = main_dialog.treeView.currentIndex()


    # use the filePath method of the file system model to get the path of the file associated with that index.
    filepath_str = main_dialog.fileSystemModel.filePath(index)

    filepath_list = os.path.normpath(filepath_str).split(os.sep)


    if "Results" not in filepath_list:
        main_dialog.log_raw_signal.emit(f"No organism was selected. Please select an organism, a subgroup, a group or a kingdom in the arborescence.")
        main_dialog.log_raw_signal.emit(f"Example: <em style=\"color: lightblue;\">Bacteria/Pseudomonadota/Gammaproteobacteria/Candidatus Carsonella ruddii</em>.")
        main_dialog.log_signal.emit("Done.", {"color":"green", "bold":True})
        main_dialog.log_raw_signal.emit("")
        return

    if selected_regions==[]:
        main_dialog.log_raw_signal.emit("No region was selected. Please select at least one functional region.")
        main_dialog.log_raw_signal.emit(f"Example: <em style=\"color: lightblue;\">intron</em>.")
        main_dialog.log_signal.emit("Done.", {"color":"green", "bold":True})
        main_dialog.log_raw_signal.emit("")
        return

    relative_filepath_as_list = filepath_list[filepath_list.index("Results")+1:]
    entity = ["kingdom", "group", "subgroup", "organism"][len(relative_filepath_as_list)-1]
    if entity in ["kingdom", "group", "subgroup"]:
        main_dialog.log_signal.emit(f"Looping through the organisms in the {entity} <em style=\"color: lightblue;\">{relative_filepath_as_list[-1]}</em>...", {"bold":True})

    #DEBUG: test si collect_leaf_repertories_path fonctionne
    #all_organisme_paths = collect_leaf_repertories_path(filepath_str)
    #if all_organisme_pathsis empty,log an error
    #if all_organisme_paths == []:
    #    main_dialog.log_raw_signal.emit("noting in all_organisme_paths")
    #    return

    #for path in all_organisme_paths:
    #    main_dialog.log_signal.emit(f'[DEBUG] Organisme path: "{path}"\n', {"color": "orange"})
    #return

    organism_paths = collect_organism_paths(filepath_str)

    #sort all_organisme_paths by alphabetical order
    organism_paths.sort()

    #transform each paths in all_organisme_paths into a list of strings with the split method
    organism_paths_as_lists = []
    for path in organism_paths:
        organism_paths_as_lists.append(os.path.normpath(path).split(os.sep))

    #filepath_list = os.path.normpath(filepath_str).split(os.sep)

    #DEBUG
    #print(filepath_list)

    # get the path from "Results" to the end for each paths list in all_organisme_paths_list
    organism_relative_paths_as_list=[]
    for path_as_list in organism_paths_as_lists:
        relative_path_as_list = path_as_list[path_as_list.index("Results")+1:]
        organism_relative_paths_as_list.append(relative_path_as_list)

        #DEBUG
        #print(relative_filepath_list)

        #pas besoin de ça je pense
        #valid_path = (len(relative_filepath_list) == 4)

        # #print the path in log
        # main_dialog.log_signal.emit(f'[DEBUG] Path to organism: "{path_to_organism}"', {"color": "orange"})

    #if "Results" not in filepath_list or not valid_path:
    #    main_dialog.log_raw_signal.emit(f"No organism was selected. Please select an organism in the arborescence. (Not a kingdom, a group or a subgroup!)")
    #    main_dialog.log_raw_signal.emit(f"Example: <em style=\"color: lightblue;\">Bacteria/Pseudomonadota/Gammaproteobacteria/Candidatus Carsonella ruddii</em>.")
    #    main_dialog.log_raw_signal.emit(f"Counterexample: <em style=\"color: lightpink;\">Bacteria/Pseudomonadota/Gammaproteobacteria</em> (this is a group, not an organism).")

    main_dialog.set_value_progress_bar1_signal.emit(0)
    main_dialog.set_max_progress_bar1_signal.emit(len(organism_relative_paths_as_list))
    main_dialog.set_format_progress_bar1_signal.emit("%v out of %m")
    progress_bar1_value = 0
    #Extract the regions specified un write it in a file, for all organisms in all_organisme_relative_paths_list
    for relative_filepath_list in organism_relative_paths_as_list:
        downloadNCRegion(relative_filepath_list, selected_regions, main_dialog)

        progress_bar1_value += 1
        main_dialog.set_value_progress_bar1_signal.emit(progress_bar1_value)
        main_dialog.set_value_progress_bar2_signal.emit(100)


        if main_dialog.processing_thread.stop_requested:
            return

    main_dialog.log_signal.emit("Done.", {"color":"green", "bold":True})
    main_dialog.log_raw_signal.emit("")

def on_stop_button_pressed(main_dialog: MainDialog):
    """This function is executed when the button Stop is pressed."""
    main_dialog.log_signal.emit("Stopped.", {"color":"orange", "bold":True})
    main_dialog.log_raw_signal.emit("")


#Recursively collect paths of all leaf files and empty directories in a directory tree.
def collect_organism_paths(root_path):
    paths = []
    for name in os.listdir(root_path):
        path = os.path.join(root_path, name)
        if os.path.isdir(path):
            child_paths = collect_organism_paths(path)
            if len(child_paths) == 0:
                paths.append(path)
            else:
                paths += child_paths
    #si paths est vide ici, ça veut dire que root path est une feuille (un organisme), donc on le met dans paths
    if paths == []:
        paths.append(root_path)
    return paths