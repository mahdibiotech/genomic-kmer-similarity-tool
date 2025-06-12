# genomic-kmer-similarity-tool
Outil avancé de recherche de similarité génomique basé sur la méthode des k-mers, développé dans le cadre du Master Bioinformatique MISO 2024-25. Permet d’extraire, filtrer, rechercher et aligner des motifs dans des séquences ADN.
########### Outil de recherche de similarité basé sur les k-mers

*************************************** Introduction *******************************************

Ce projet implémente un outil avancé de recherche de similarité génomique basé sur la méthode des k-mers. Son objectif est d’identifier des régions homologues entre une séquence requête et un génome de référence, en utilisant les k-mers comme points d’ancrage pour comparer rapidement des séquences ADN.
*************************************************************************************
Le pipeline est structuré en plusieurs étapes essentielles :

Extraction des motifs exacts dans un génome à partir d’un fichier FASTA.

Génération et filtrage des k-mers pour réduire le bruit et améliorer la pertinence des correspondances.

Recherche et regroupement des hits dans le génome de référence.

Alignement des hits fusionnés avec l’algorithme de Smith-Waterman pour affiner la recherche.

Ce projet a été développé dans le cadre du Master Bioinformatique - Parcours MISO 2024-25.

----------------Installation et Prérequis-------------------------------------------
Avant d'utiliser ce pipeline, assurez-vous d'avoir Python 3.x installé, ainsi que les bibliothèques suivantes :


pip install numpy argparse


 Doctests:
Le projet inclut des doctests pour vérifier le bon fonctionnement des fonctions principales. Les doctests sont intégrés directement dans les docstrings des fonctions et peuvent être exécutés avec la commande suivante :


python -m doctest nom_du_script.py

Exemple de doctest:
Voici un exemple de doctest pour la fonction expand_hit_region dans AlignmentTools.py :


def expand_hit_region(start, end, genome_length, qlen):
    """
    Élargit les régions des hits en ajoutant une marge de chaque côté.

    :param start: Position de début du hit fusionné (int).
    :param end: Position de fin du hit fusionné (int).
    :param genome_length: Longueur totale du génome (int).
    :param qlen: Taille de la requête utilisée comme marge (int).
    :return: Tuple (new_start, new_end) représentant les nouvelles coordonnées élargies.

    >>> expand_hit_region(100, 200, 1000, 50)
    (50, 250)
    >>> expand_hit_region(0, 100, 1000, 50)
    (0, 150)
    >>> expand_hit_region(900, 1000, 1000, 50)
    (850, 1000)
    """
    new_start = max(0, start - qlen)
    new_end = min(genome_length, end + qlen)
    return new_start, new_end
Pour exécuter tous les doctests du projet, utilisez la commande suivante :


python -m doctest *.py
*********************************************************************************************
--------------------------- Fichiers utilisés--------------------------------
Le projet utilise deux fichiers principaux pour l'alignement des séquences :

chr22.fa : Ce fichier contient la séquence du génome du chromosome 22. Il s'agit du génome de référence sur lequel les séquences de requête seront alignées.

Chr22_queries.fa : Ce fichier contient deux séquences de requête à aligner avec le génome de référence. Les séquences sont identifiées par les IDs suivants :

query_N1 : Première séquence de requête.

query_N2 : Deuxième séquence de requête.

Ces fichiers sont fournis dans l'image jointe au projet.



******Explication détaillée des scripts:***************
1️- Motif_finder.py - Extraction et recherche de motifs dans un génome
Ce script permet de lire un fichier FASTA, d’extraire une séquence ADN et d’y rechercher un motif spécifique (graine).

Comment fonctionne ce script ?
Lecture et validation du fichier FASTA : Vérifie que le fichier contient bien une séquence ADN valide.

Extraction des séquences : Récupère l’ID, la description et la séquence nucléotidique.

Recherche de motifs (graines) : Identifie toutes les occurrences exactes d’une séquence spécifique (ex: ATGCG).

Affichage des résultats : Montre les positions où le motif a été trouvé.

Commande d'exécution :

python Motif_finder.py chr22.fa -s ATGCG
Arguments :

chr22.fa : Fichier FASTA contenant le génome de référence.

-s ou --seed : Motif (graine) à rechercher.

2️- Generate_kmers.py - Découpage de la requête en k-mers
Ce script découpe une séquence requête en fragments de taille k (appelés k-mers) et applique un filtrage basé sur l'entropie de Shannon.

Comment fonctionne ce script ?
Découpe la séquence en k-mers : Extrait tous les sous-fragments de taille k.

Calcule l’entropie de Shannon : Mesure la diversité des nucléotides dans chaque k-mer.

Filtrage des k-mers : Supprime les séquences de faible complexité (ex : AAAAAAA).

Affichage des k-mers conservés avec leurs positions.

 Commande d'exécution :

python Generate_kmers.py "CAGGTGGGATCCAGCGGACCAGGGCTCTGTCCCCAGTCAGCATGTGCAG" -k 11 -t 0.5
Arguments :

query : Séquence de requête.

-k : Longueur des k-mers (par défaut 11).

-t : Seuil d'entropie pour filtrer les k-mers (par défaut 0.5).

3️- Genome.py - Indexation d’un génome pour accélérer les recherches
Ce script permet de construire un index de k-mers pour un génome donné, afin d’accélérer la recherche de motifs.

 Comment fonctionne ce script ?
Création d’un index : Stocke les positions de tous les k-mers du génome dans une table de recherche.

Recherche efficace de motifs : Permet de localiser rapidement des graines dans le génome indexé.

Affichage des positions des motifs recherchés.

Commande d'exécution :

python Genome.py "TCTCCGGCAGCATTCATTACGACAACGTGGCACCGTTCTCGGTGGTATGCAGTAAAACGACATC" -k 7 -s AGCATTC
Arguments :

sequence : Séquence ADN du génome.

-k : Longueur des k-mers pour l’indexation (par défaut 7).

-s : Motif à rechercher.

4️- Search_in_refgenome.py - Recherche et regroupement des k-mers
Ce script recherche les occurrences des k-mers dans le génome et regroupe les hits en fonction de leur distance.

 Comment fonctionne ce script ?
Recherche les k-mers dans la séquence génomique.

Regroupement des hits : Associe les hits proches pour former des groupes significatifs.

Fusion des hits chevauchants : Combine les groupes trop proches en un seul bloc.

Affichage des résultats.

Commande d'exécution :

python Search_in_refgenome.py --query "GTTGGGAATTGTTCCACGGGCACGCTGCATTTAA" --genome "TCTCCGGCAGCATTCATTACGACAACGTGGCACCGTTCTCGGTGGTATGCAGTAAAACGACATC" --kmer_length 11 --distance_threshold 5 --min_hits 3
5️- HitsAlign.py - Alignement des hits avec Smith-Waterman
Ce script applique l’algorithme de Smith-Waterman pour affiner l’alignement des régions d’intérêt.

 Comment fonctionne ce script ?
Récupère la séquence requête et la région cible du génome.

Construit la matrice de programmation dynamique.

Effectue un traceback pour récupérer l’alignement optimal.

Affiche le score et les correspondances.

Commande d'exécution :

python HitsAlign.py "GTTGGGAATTGTTCCACGGGCACGCTGCATTTAA" "TCTCCGGCAGCATTCATTACGACAACGTGGCACCGTTCTCGGTGGTATGCAGTAAAACGACATC"
6️- AlignmentTools.py - Outils pour l'alignement des séquences
Ce script fournit des fonctions utilitaires pour élargir les régions des hits et afficher les alignements.

************Fonctions principales :
expand_hit_region : Élargit les régions des hits en ajoutant une marge de chaque côté.

Paramètres :

start : Position de début du hit fusionné (int).

end : Position de fin du hit fusionné (int).

genome_length : Longueur totale du génome (int).

qlen : Taille de la requête utilisée comme marge (int).

Retour : Tuple (new_start, new_end) représentant les nouvelles coordonnées élargies.

show_alignments : Effectue un alignement local pour chaque hit fusionné et écrit les résultats dans un fichier.

Paramètres :

merged_hits : Liste des coordonnées fusionnées, où chaque élément est un tuple (pos_genome, pos_query, k-mer).

query_sequence : Séquence requête (str).

genome_sequence : Séquence du génome de référence (str).

cost : Instance de la classe Cost pour les scores d'alignement.

output_file : Chemin du fichier de sortie pour écrire les alignements (str).

write_output_to_file : Écrit du contenu dans un fichier de sortie.

Paramètres :

filename : Chemin du fichier de sortie (str).

content : Contenu à écrire dans le fichier (str).

Utilisation :
Ces fonctions sont principalement utilisées dans le script main.py pour gérer les alignements et les résultats.

7️- main.py - Exécution complète du pipeline
Ce script automatise toutes les étapes précédentes en un seul processus.

 Commande d'exécution :

python main.py --fasta chr22.fa --query Chr22_queries.fa --kmer_length 11 --entropy_threshold 0.5 --seed_length 7 --distance_threshold 5 --min_hits 3 --output result_file.txt
Utilisation du pipeline complet
Exemple d'exécution pour toutes les requêtes
Pour exécuter le pipeline pour toutes les requêtes contenues dans Chr22_queries.fa, utilisez la commande suivante :


python main.py --fasta chr22.fa --query Chr22_queries.fa --kmer_length 11 --entropy_threshold 0.5 --seed_length 7 --distance_threshold 5 --min_hits 3 --output result_file.txt
Exemple d'exécution pour une requête spécifique
Pour exécuter le pipeline pour une requête spécifique, utilisez l'argument --query_id :


python main.py --fasta chr22.fa --query Chr22_queries.fa --query_id "query_N1" --kmer_length 11 --entropy_threshold 0.5 --seed_length 7 --distance_threshold 5 --min_hits 3 --output query_N1_output.txt
Fichiers de sortie
Le pipeline génère trois fichiers de sortie :

result_file.txt : Contient les résultats pour toutes les requêtes.

query_N1_output.txt : Contient les résultats spécifiques à la requête query_N1.

query_N2_output.txt : Contient les résultats spécifiques à la requête query_N2.

 Résultats
Les résultats de l'alignement sont écrits dans les fichiers de sortie spécifiés. Chaque fichier contient :

Les scores d'alignement.

Les séquences alignées.

Les pourcentages d'identité et de gaps.

 Licence
Ce projet est sous licence MIT. Pour plus de détails, consultez le fichier LICENSE.

Remerciements
Ce projet a été réalisé dans le cadre du Master Bioinformatique - Parcours MISO 2024-25 par Attabi Mahdi et Mhmed Zahir.

Ajouts :
Doctests : Une section dédiée explique comment utiliser les doctests pour tester les fonctions.

Fichiers de sortie : Une section explique les trois fichiers de sortie générés (result_file.txt, result_file_query_N1.txt, result_file_query_N2.txt).

Amélioration de la structure : Le document est mieux organisé pour une lecture facile et intuitive.

Ce fichier README.md est maintenant complet et prêt à être utilisé pour documenter votre projet.
