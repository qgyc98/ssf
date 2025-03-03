#!/bin/bash
# created on 14th of August 2008 in EPCC by JB
# broz@cml.fsv.cvut.cz
# pro beh programu je nezbytna instalace balicku zenity
# skript provede rozrezani site NAME.top pomoci programu metisprep, partdmesh a metispost na casti

# pocatecni inicializace globalnich promnenych
ZENITY=""
odpoved=""
NAME_SIT=""
ND=""
MD=""
DIR=""
KOPIE=""
NAME_OUT=""

# pocet parametru
poc_param=$#
# prepis jednotlivych parametru
par1=$1
par2=$2
par3=$3
par4=$4
# funkce provede kontrolu zda je nainstalovany balicek zenity, ktery je pouzit pro GUI. Predpoklada se umisteni v /bin/bash
kontrola_zenity (){
    echo "Pro spravny beh skriptu je nutne mit balicek zenity"
    sleep 1
    if [  -e /usr/bin/zenity ]
    then
	ZENITY=1
	echo "Balicek zenity je naistalovan v /usr/bin/"
    else
	ZENITY=0
	echo "Balicek zenity neni nainstalovan v /usr/bin/ bude nutne zadavat vse z radky"
    fi
}

# funkce zkontroluje zda je zadan prvni parametr, pokud ne je zahajeno interaktivni spusteni
kontrola_spusteni(){
    if  [ $poc_param -eq 0 ]
    then 
	if [ $ZENITY -eq 1 ]
	then
	    odpoved=$(zenity --question --title="Interaktivni spusteni" --text "Pro spusteni skriptu $0 jste zadali malo parametru, chcete zahajit interaktivni spusteni?"; echo $?)
	    case "$odpoved" in
		0 )
		    echo "Zahajeno interaktivni spusteni"
		    ;;
		1 )
		    echo "Skript bude ukoncen"
		    exit 1
		    ;;
		* )
		    echo "Spatna volba, skript bude ukoncen"
		    exit 1
		    ;;
	    esac
	else
	    echo  "Malo parametru, chcete zahajit interaktivni spusteni?"
	    echo "Odpoved:Ano, Ne"
	    read odpoved
	    case "$odpoved" in
		[aA] | [aA][nN][oO] )
		    echo "Zahajeno interaktivni spusteni"
		    ;;
		[nN] | [nN][eE] )
		    echo "Skript bude ukoncen"
		    echo "Navod na spusteni:$0 jmeno_pro_spusteni pocet_domen_v_x pocet_domen_v_y sila_v_x sila_v_y"
		    echo "Priklad:$0 boundparcg 2 2 100000 0"
		exit 1
		;;
		* )
		    echo "Spatna volba, skript bude ukoncen"
		    echo "Navod na spusteni:$0 jmeno_pro_spusteni pocet_domen_v_x pocet_domen_v_y sila_v_x sila_v_y"
		    echo "Priklad:$0 boundparcg 2 2 100000 0"
		    exit 1
		    ;;
	    esac
	fi
    else
	echo "Neinteraktivni spusteni"
    fi
}

# funkce doplni jmeno site pro rozrezani
doplneni_jmena_site (){
# doplneni jmena souboru
    if [ $poc_param -eq 0 ]
    then 
	if [ $ZENITY -eq 1 ]
	then
	    NAME_SIT=$(zenity --file-selection --title="Vyber souboru - Vyberte prosim soubor se siti pro rozrezani")
	else
	    echo "Zadejte, prosim, jmeno souboru, ktery ma byt rozrezan" 
	    read NAME_SIT
	fi
    else
	NAME_SIT=$par1
    fi
}

# funkce zjisti zda bylo nejake jmeno zadano, pokud ne, zavola funkci na doplneni
kontrola_jmena_site (){
# test zda je opravdu zadan soubor se siti
    if [ -z $NAME_SIT ]
    then
	if [ $ZENITY -eq 1 ]
	then
 	    zenity --error --text "Nebyl zadan soubor se siti pro rozrezani, je nutne jmeno zadat"
	    doplneni_jmena_site
	    if [ -z $NAME_SIT ]
	    then
		zenity --error --text "Opet nebyl zadan soubor se siti pro rozrezani,skript bude ukoncen"
		exit 1
	    fi
	else
	    echo "Nebyl zadan soubor se siti pro rozrezani"
	    echo "Jmeno souboru se siti je nutne zadat"
	    doplneni_jmena_site
	    if [ -z $NAME_SIT ]
	    then
		zenity --error --text "Opet nebyl zadan soubor se siti pro rozrezani,skript bude ukoncen"
		exit 1
	    fi
	fi
    else
	echo "Byl zadan nasledujici soubor $NAME_SIT"
    fi 
}

# funkce zkontroluje existenci souboru a jeho citelnost
test_souboru_site () {
#test existence souboru se siti
    echo "Za 1 sekundu bude spusten test existence souboru $NAME_SIT"
    if [ -e $NAME_SIT ]
    then
	echo "Soubor $NAME_SIT existuje"
    #testovani citelnosti souboru
	if [  -r $NAME_SIT ]
	then 
	    echo "Soubor $NAME_SIT je citelny"
	else
	    echo "Soubor $NAME_SIT neni citelny"
	fi
    else
	zenity --error --text "Soubor $NAME_SIT neexistuje"
	NAME_SIT=$(zenity --file-selection --title="Vyber souboru - Vyberte prosim soubor se siti pro rozrezani")
    # test existence a citelnosti dodatecneho souboru se sablonou
	if [ -e $NAME_SIT ]
	then
	    echo "Soubor $NAME_SIT existuje"
    #testovani citelnosti souboru
	    if [  -r $NAME_SIT ]
	    then 
		echo "Soubor $NAME_SIT je citelny"
	    else
		echo "Soubor $NAME_SIT neni citelny"
	    fi
	else
	    NAME_SIT=$(zenity --file-selection --title="Vyber souboru" --text "Vyberte prosim soubor se siti pro rozrezani")
	fi
    fi
}

# funkce doplni pocet domen pro rozrezani
doplneni_poctu_domen () {
# doplneni poctu domen
    if [  $poc_param  -lt 2 ]
    then
	if [ $ZENITY -eq 1 ]
	then
	    ND=$(zenity --entry --title="Pocet domen" --text "Zadejte prosim pocet domen na ktere chcete sit rozrezat")
	else
	    echo "Zadejte prosim pocet domen na ktere chcete sit rozrezat"
	    read ND
	fi
    else
	ND=$par2
    fi
}

#funkce otestuje zadani poctu domen
test_poctu_domen () {
#test zda byl skutecne zadan pocet domen
    if [ -z $ND ]
    then
	if [ $ZENITY -eq 1 ]
	then
 	    zenity --error --text "Nebyl zadan pocet domen, pocet je nutne zadat"
	    doplneni_poctu_domen
	    if [ -z $ND ]
	    then
 		zenity --error --text "Opet nebyl zadan pocet domen, skript bude ukoncen"
		exit 1
	    fi
	else
	    echo "Nebyl zadan pocet domen pro rozrezani, pocet je nutne zadat"
	    doplneni_ND
	    if [ -z $ND ]
	    then
		echo "Opet nebyl zadan pocet domen, skript bude ukoncen"
		exit 1
	    fi
	fi
    else
	echo "Byl zadan tento pocet domen: $ND"
    fi
}

#funkce doplni mesh description
doplneni_MD () {
# doplneni Mesh Description
    if [  $poc_param  -lt 3 ]
    then
	if [ $ZENITY -eq 1 ]
	then
	    MD=$(zenity --list  --text "Vyberte, prosim, mesh description" --radiolist  --column "Vyber" --column "Mesh description" TRUE "All nodes"  FALSE "Boundary nodes" )
	else
	    echo "Vyberte, prosim, mesh description"
	    echo "1 - All nodes" 
	    echo "2 - Boundary nodes"
	    read MD
	fi
	case "$MD" in
	    "All nodes" )
		MD=1
		echo "Jako mesh description bude pouzita volba all nodes"
		;;
	    "Bound nodes" )
		MD=2
		echo "Jako mesh description bude pouzita volba bound nodes"
		;;
	esac
    else
	MD=$par3
    fi
}

#funkce otestuje zda bylo zadano mesh description
test_MD () {
#test zda byl opravdu zadan mesh description
    if [ -z $MD ]
    then
	if [ $ZENITY -eq 1 ]
	then
	    zenity --error --text "Nebyl zadan mesh description, skript bude ukoncen"
 	    doplneni_MD
	    if [-z $MD ]
	    then
		zenity --error --text "Opet nebyl zadan mesh description, skript bude ukoncen"
		exit 1
	    fi
	else
	    echo "Nebyl zadan mesh description, skript bude ukoncen"
 	    doplneni_MD
	    if [-z $MD ]
	    then
		echo "Opet nebyl zadan mesh description, skript bude ukoncen"
		exit 1
	    fi
	fi
    else
	echo "Mesh description zadan"
    fi
}

#funkce zjisti zda je v aktualnim adresari program metisprep, pokud ne spusti funkci na vyber adresare a jeho zkopirovani
test_metisprep () {
#test existence metisprep
    echo "Za 1 sekundu bude spusten test programu metisprep"
    if [ -e metisprep ]
    then
	echo "Program metisprep je v aktualnim adresari $PWD"
    #testovani spustitelnosti
	if [  -x metisprep ]
	then 
	    echo "Program metisprep je spustitelny"
	else
	    echo "Program metisprep neni spustitelny"
	    chmod a+x metisprep
	fi
    else
	zenity --error --text "Program metisprep neni v adresari $PWD, predpokladene umisteni SIFEL/PREP/PARTITIONING"
	#kopirovani_metisprep
	if [ -z $DIR ]
	then
	    DIR=$(zenity --file-selection --directory --title="Vyber umisteni - Vyberte prosim umisteni programu metisprep")
	fi
	cp $DIR/metisprep .
	echo "Program metisprep byl nakopirovan z adresare $DIR do adresare $PWD"
	KOPIE=metisprep
# test existence
	if [ -e metisprep ]
	then
	    echo "Program metisprep je v adresari $PWD"
	        #testovani spustitelnosti
	    if [  -x metisprep ]
	    then 
		echo "Program metisprep je spustitelny"
	    else
		echo "Program metisprep neni spustitelny"
		chmod a+x metisprep
	    fi
	fi
	echo "****************************************************************************"
    fi
}

#funkce zjisti zda je v aktualnim adresari program partdmesh, pokud ne spusti funkci na vyber adresare a jeho zkopirovani
test_partdmesh () {
#test existence partdmesh
    echo "Za 1 sekundu bude spusten test programu partdmesh"
    if [ -e partdmesh ]
    then
	echo "Program partdmesh je v adresari $PWD"
#testovani citelnosti souboru
	if [  -x partdmesh ]
	then 
	    echo "Program partdmesh je spustitelny"
	else
	    echo "Program partdmesh neni spustitelny"
	    chmod a+x partdmesh
	fi
    else
	if [ -z $DIR ]
	then
	    zenity --error --text "Program partdmesh neni v adresari $PWD, predpokladene umisteni SIFEL/PREP/PARTITIONING"
	    DIR=$(zenity --file-selection --directory --title="Vyber umisteni - Vyberte prosim umisteni programu partdmesh")
	    KOPIE=partdmesh
	else
	    zenity --error --text "Program partdmesh neni v adresari $PWD a bude nakopirovan z adresare $DIM, stejne jako program $KOPIE"
	    KOPIE=metisprep
	fi
	cp $DIR/partdmesh .
	
	echo "Program partdmesh byl nakopirovan z adresare $DIR do adresare $PWD"
# test existence a citelnosti dodatecneho souboru se sablonou
	if [ -e partdmesh ]
	then
	    echo -n "Program partdmesh je v adresari $PWD"
	    pwd
    #testovani citelnosti souboru
	    if [  -x partdmesh ]
	    then 
		echo "Program partdmesh je spustitelny"
	    else
		echo "Program partdmesh neni spustitelny"
		chmod a+x partdmesh
	    fi
	fi
	echo "****************************************************************************"
    fi
}

#funkce zjisti zda je v aktualnim adresari program metispost, pokud ne spusti funkci na vyber adresare a jeho zkopirovani
test_metispost () {
#test existence metispost
    echo "Za 1 sekundu bude spusten test programu metispost"
    if [ -e metispost ]
    then
	echo "Program metispost je v adresari $PWD"
	pwd
    #testovani citelnosti souboru
	if [  -x metispost ]
	then 
	    echo "Program metispost je spustitelny"
	else
	    echo "Program metispost neni spustitelny"
	    chmod a+x metispost
	fi
    else
	if [ -z $DIR ]
	then
	    zenity --error --text "Program metispost neni v adresari $PWD, predpokladene umisteni SIFEL/PREP/PARTITIONING"
	    DIR=$(zenity --file-selection --directory --title="Vyber umisteni - Vyberte prosim umisteni programu metispost")
	else
	    zenity --error --text "Program metispost neni v adresari $PWD a bude nakopirovan z adresare $DIR stejne jako $KOPIE"
	fi
	cp $DIR/metispost .
	echo "Program metispost byl nakopirovan z adresare $DIR do adresare $PWD"
# test existence a citelnosti dodatecneho souboru se sablonou
	if [ -e metispost ]
	then
	    echo  "Program metispost je v aktualnim adresari $PWD"
	    #testovani citelnosti souboru
	    if [  -x metispost ]
	    then 
		echo "Program metispost je spustitelny"
	    else
		echo "Program metispost neni spustitelny"
		chmod a+x metispost
	    fi
	fi
	echo "****************************************************************************"
    fi
}

#funkce doplni vystupni jmeno pro rozrezane soubory
doplneni_vystupniho_jmena () {
# doplneni vstupniho jmena
    if [  $poc_param  -lt 4 ]
    then
	if [ $ZENITY -eq 1 ]
	then
	    NAME_OUT=$(zenity --entry --title="Jmeno vystupnich soubour" --text "Zadejte, prosim, jmeno vystupnich souboru, bez pripony")
	else
	    echo "Zadejte, prosim, jmeno vystupnich souboru, bez pripony"
	    read NAME_OUT
	fi
    else
	NAME_OUT=$par4
    fi
}

#funkce otestuje zda bylo zadano nejake vystupni jmeno
test_vystupniho_jmena () {
#test zadani vystupniho jmena
    if [ -z $NAME_OUT ]
    then
	if [ $ZENITY -eq 1 ]
	then
 	zenity --error --text "Nebylo zadano jmeno vystupnich souboru, jmeno je nutne zadat"
	doplneni_vystupniho_jmena
	if [ -z $NAME_OUT ]
	then
 	    zenity --error --text "Opet nebylo zadano jmeno vystupnich souboru, skript bude ukoncen"
	    exit 1
	fi
	else
	    echo "Nebylo zadano jmeno vystupnich souboru, jmeno je nutne zadat"
	    doplneni_vstupniho_jmena
	    if [ -z $NAME_OUT ]
	    then
 		echo "Opet nebylo zadano jmeno vystupnich souboru, skript bude ukoncen"
		exit 1
	    fi
	fi
    else
	echo "Zadano jmeno $NAME_OUT pro vystupni soubory"
    fi
}

#hlavni cast skriptu
kontrola_zenity
kontrola_spusteni
doplneni_jmena_site
kontrola_jmena_site
test_souboru_site
echo "****************************************************************************"
echo "Jmeno pro vypocet je $NAME_SIT"
sleep 1
echo "****************************************************************************"
doplneni_poctu_domen
test_poctu_domen
doplneni_MD
test_MD
test_metisprep
echo "========================================================"
echo "PREPROCESING PRO PARTDMESH"
./metisprep $NAME_SIT tmp.txt
echo "========================================================"
echo "METIS"
echo "========================================================"
test_partdmesh
./partdmesh tmp.txt $ND
echo "========================================================"
echo "POSTPROCESING"
test_metispost
doplneni_vystupniho_jmena
test_vystupniho_jmena
./metispost $NAME_SIT tmp.txt.epart.$ND $ND $NAME_OUT $MD
echo "========================================================"
# uklid
rm tmp.* 
 