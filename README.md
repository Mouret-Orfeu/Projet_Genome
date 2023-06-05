# Projet Genbank

## Principe du projet
Afin de récupérer les séquences d'ADN d'une espèce donnée de la base de données Genbank, nous avons développé une interface graphique Python affichant l'arbre du vivant des espèces disponibles sous forme d'arborescence de dossiers. 
L'interface permet de selectionner les espèces en cliquant sur le dossier correspondant.
Ensuite il suffit de lancer la recherche et le parsing des données de Genbank avec le bouton start, ce qui remplira les dossiers choisis avec les fichiers contenants les types de séquences d'ADN cochées au préalables (CDS et/ou intron).

Ce projet a été réalisé par :
- Maxime Drouhin
- Maty Laye Souare
- Orfeu Mouret

## Prérequis

- Python 3.8+ (versions antérieures non testées)

## Exécution

Double-cliquer sur `run.bat` pour lancer le programme. Cela installe les packages requis et exécute le programme.

### Exécution manuelle

- Installer les packages requis avec la commande `pip install -r requirements.txt`.

Puis il suffit d'exécuter `genome-project.py` avec Python avec un simple double-clic ou avec la commande suivante :

Sous Windows :
```
python genome-project.py
```

Sous Linux :
```
python3 genome-project.py
```

## Nettoyer les fichiers générés

Double-cliquer sur `clean.bat` pour supprimer les fichiers générés. Cela supprime les dossiers `data` et `Results`.

### Nettoyage manuel

Il suffit de supprimer les dossiers `data` et `Results`.

## Informations concernant les fonctionnalités

- On ne peut pour l'instant que sélectionner des organismes, sous-groupes, groupes et règnes. Pas de multi-sélection possible.

## Implémentation

- Nous avons choisi, pour la robustesse du code, que le bouton "Stop" ne coupe pas immédiatement le thread de processing. Au lieu de ça, appuyer sur "Stop" met un flag, et le thread de processing surveille ce flag. Il peut donc s'arrêter proprement quand l'utilisateur demande l'arrêt, après avoir fini le parsing de la région en cours de traitement. La limite de cette approche est qu'il faut parfois attendre un peu pour que ça s'arrête.
- `overview.txt` (https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/overview.txt) est téléchargé automatiquement au premier lancement du programme. Pour le mettre à jour, ne pas hésiter à suivre les instructions de la section `Nettoyer les fichiers générés` pour qu'il soit retéléchargé.
- Pour assurer la robustesse du code, nous avons utilisé le typage.
  ```py
  def downloadNCRegion(organism_path_list: List[str], regions: str, main_dialog: MainDialog):
    # ...
  ```
- Pour ouvrir les fichiers de manière sécurisée, nous avons utilisé des gestionnaires de contexte. Nous utilisons aussi le standard `utf-8` pour une meilleure compatibilité.
  ```py
  with open(file_path, "r", encoding='utf-8') as f:
      # ...
  ```
  ```
- Compatibilité avec divers systèmes d'exploitation grâce à l'utilisation de `os.path`.
  Par exemple, `os.path.join` permet de créer des chemins de fichiers compatibles avec les systèmes d'exploitation Windows et Linux, et évite d'utiliser la syntaxe `path/to/file` qui peut ne pas fonctionner sur Windows.
  ```py
  path = os.path.join("Results", kingdom, group, subgroup, organism)
  ```
