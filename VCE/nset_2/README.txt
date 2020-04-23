Makefile ci joint :
make pour compiler
make run pour compiler (si necessaire) et executer
make time, pareil que run mais avec la mesure de temps avec time ./exec en plus
make clean pour delete l'executable

3 arguments a preciser, mais avec une valeur par defaut si non precisee
file : chemin de la matrice a utiliser (les donnees)
size : valeur de m a utiliser (12 par defaut)
nbr_threads : nombre de threads openmp a utiliser (18 par defaut, c'est le max sur un seul noeud sur gpm)
