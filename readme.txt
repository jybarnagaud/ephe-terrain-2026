--- EPHE - master Biodiversité Environnement --- 

UE terrain : dossier du projet écoacoustique
sept 2026
encadrement : Jean-yves.barnagaud@ephe.psl.eu
étudiants
Sarah Merle
Lou Sausse Corbère
Killian Rey
Mila Lefèvre

!!! les enregistrements audio d'origine ne sont pas sauvegardés dans 
ce dossier de projet (9Go). sauvegardes en local sur support physique, 
à demander à l'encadrant
------

# 1 objectif
------------

étude exploratoire de la variation du paysage acoustique sur un gradient de paysages ouverts - forestiers en milieu agropastoral de moyenne montagne

# 2 contexte
------------

lieu : nord-Vercors (Isère, Fr)
date : septembre 2026
matériel : Wildlife Acoustics SM Mini (1ère génération) avec 1 micro omnidirectionnel dans l'audible
software : toutes les analyses sont réalisées sous R. Acquisition des points GPS par Organic Maps sous Android (smartphone)

# 3 contenu du dossier
----------------------

> racine
readme : ce fichier
ephe-terrain-2026 : batchfile projet RStudio
.Rhistory
.gitignore

>> anciens-rapports : rapports des années précédentes
>> bibliographie : références majeures
>> data
	>>> sig : coordonnées WGS84 des enregistreurs + métadonnées
>> outputs : toutes les sorties des scripts R
>> r-functions : fonctions R pour le calcul des indices écoacoustiques
>> scripts
		>>>> ephe-terrain-2026-indices-ecoacoustiques.R : calcul des indices écoacoustiques à partir des fichiers audio (sous R)
		>>>> ephe-terrain-2026-analyses.R : analyse des données (sous R - nécessite run des indices écoacoustiques en amont)

