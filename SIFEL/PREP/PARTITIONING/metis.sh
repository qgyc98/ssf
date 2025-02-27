#!/bin/bash
# JB 190407
# skript provede rozrezani site NAME.top pomoci programu metisprep, partdmesh a metispost na casti
echo "Zadejte prosim nazev bez pripon souboru ktery chcete rozrezat - musi mit priponu .top"
read NAME
# testovani zda soubor opravdu existuje
if test -f $NAME.top
# Soubor existuje
then echo "Bylo zadano toto jmeno: $NAME.top"
else
# Soubor neexistuje.
echo "Spatne zadane jmeno souboru, zadejte jej prosim znovu"
read NAME   
fi
echo "Zadejte prosim pocet casti na ktere chcete rozrezat sit"
read ND
echo "Bylo zadan tento pocet: $ND"
echo "Zadejte prosim mesh description pro rozrezanou sit"
echo "1 - all nodes"
echo "2 - bound nodes"
read MD
case $MD in
  1)
  echo "Jako mesh description bude pouzita volba all nodes";;
  2)
  echo "Jako mesh description bude pouzita volba bound nodes";;
esac

echo" ========================================================"
echo "PREPROCESING PRO PARTDMESH"
#echo ./metisprep $NAME.top $NAME.txt
 ./metisprep $NAME.top $NAME.txt
echo "========================================================"
echo "METIS"
echo "========================================================"
#echo ./partdmesh $NAME.txt $ND
./partdmesh $NAME.txt $ND
echo "========================================================"
echo "POSTPROCESING"
#echo ./metispost $NAME.top $NAME.txt.epart.$ND $ND $NAME $MD
./metispost $NAME.top $NAME.txt.epart.$ND $ND $NAME $MD
echo "========================================================"
 
