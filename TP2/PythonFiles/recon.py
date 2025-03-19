#!/usr/bin/env python
# -*- coding: utf-8 -*-
# TP reconstruction TDM (CT)
# Prof: Philippe Després
# programme: Dmitri Matenine (dmitri.matenine.1@ulaval.ca)


# libs
import numpy as np
import time

# local files
import geo as geo
import util as util
import CTfilter as CTfilter

## créer l'ensemble de données d'entrée à partir des fichiers
def readInput():
    # lire les angles
    [nbprj, angles] = util.readAngles(geo.dataDir+geo.anglesFile)

    print("nbprj:",nbprj)
    print("angles min and max (rad):")
    print("["+str(np.min(angles))+", "+str(np.max(angles))+"]")

    # lire le sinogramme
    [nbprj2, nbpix2, sinogram] = util.readSinogram(geo.dataDir+geo.sinogramFile)

    if nbprj != nbprj2:
        print("angles file and sinogram file conflict, aborting!")
        exit(0)

    if geo.nbpix != nbpix2:
        print("geo description and sinogram file conflict, aborting!")
        exit(0)

    return [nbprj, angles, sinogram]


## reconstruire une image TDM en mode rétroprojection
def laminogram():
    
    [nbprj, angles, sinogram] = readInput()

    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))

    center_voxel = (geo.nbvox - 1) / 2
    center_pixel = (geo.nbpix - 1) / 2
    scale = geo.voxsize / geo.pixsize
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    for j in range(geo.nbvox): # colonnes de l'image
        print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
            for a in range(len(angles)):
                x = j - center_voxel
                y = i - center_voxel
                u = -x*np.cos(angles[a]) + y*np.sin(angles[a])
                u *= scale
    
                pixel_index = round(u + center_pixel)

                image[i, j] += sinogram[a, pixel_index]
            
                
                
                
                
    util.saveImage(image, "lam")


## reconstruire une image TDM en mode retroprojection filtrée
def backproject():
    
    [nbprj, angles, sinogram] = readInput()
    
    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    
    ### option filtrer ###
    CTfilter.filterSinogram(sinogram)
    ######
    
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    for j in range(geo.nbvox): # colonnes de l'image
        print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
            for a in range(len(angles)):
                theta = angles[a]
                x = (i - geo.nbvox/2) * geo.voxsize
                y = (j - geo.nbvox/2) * geo.voxsize
                # Calculer la position du voxel projeté sur le détecteur
                s = x * np.cos(theta) + y * np.sin(theta)
                # Convertir la position en indice de détecteur
                detector_index = int((s / geo.detector_width + 0.5) * geo.nbdet)
                # Ajouter la contribution du sinogramme filtré à l'image
                if 0 <= detector_index < geo.nbdet:
                    image[i, j] += sinogram[a, detector_index]
    
                #votre code ici
               #pas mal la même chose que prédédemment
            #mais avec un sinogramme qui aura été préalablement filtré
    
    util.saveImage(image, "fbp")


## reconstruire une image TDM en mode retroprojection
def reconFourierSlice():
    
    [nbprj, angles, sinogram] = readInput()

    # initialiser une image reconstruite, complexe
    # pour qu'elle puisse contenir sa version FFT d'abord
    IMAGE = np.zeros((geo.nbvox, geo.nbvox), 'complex')
    
    # conteneur pour la FFT du sinogramme
    SINOGRAM = np.zeros(sinogram.shape, 'complex')

    #image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    #votre code ici
   #ici le défi est de remplir l'IMAGE avec des TF des projections (1D)
   #au bon angle.
   #La grille de recon est cartésienne mais le remplissage est cylindrique,
   #ce qui fait qu'il y aura un bon échantillonnage de IMAGE
   #au centre et moins bon en périphérie. Un TF inverse de IMAGE vous
   #donnera l'image recherchée.

   
    
    util.saveImage(image, "fft")


## main ##
start_time = time.time()
laminogram()
#backproject()
#reconFourierSlice()
print("--- %s seconds ---" % (time.time() - start_time))

