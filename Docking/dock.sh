#!/usr/bin/env dash

set -e
tabs -4
IFS="
"

exitWithMsg() {
	echo "$1"
	exit 1
}

#Zmienne środowiskowe
	#Lokalizacja tego skryptu
ROOTDIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
	#Ścieżka do viny
VINA=$ROOTDIR/vina
	#Ścieżka do receptora
RECEPTOR=$ROOTDIR/r.pdbqt
	#Ścieżka do folderu z ligandami
LIGDIR=$ROOTDIR/lig
	#Nazwa plików PDBQT ligandów
LIGFILE=l.pdbqt

[ -x "$VINA" ] || exitWithMsg "Nie znaleziono viny."
[ -f "$RECEPTOR" ] || exitWithMsg "Nie znaleziono pliku PDBQT receptora."
[ -d "$LIGDIR" ] || exitWithMsg "Nie znaleziono folderu z ligandami."


#Domyślne parametry
e=32			#Exhaustiveness
n=30			#Number of modes
er=100			#Energy ratio
cpu=4			#Liczba rdzeni
size=20			#Rozmiar przestrzeni dokowania


usage () {
	echo "Sposób użycia: $0 [opcje]

OPCJE															<domyślne wartości>
-e	N		Parametr 'exhaustiveness' (E)						<$e>
-n	N		Parametr 'number of modes' (N)						<$n>
-r	N		Parametr 'energy ratio'	(ER)						<$er>
-s	N		Rozmiar przestrzeni dokowania w angstremach (S)		<$size>
-c	N		Z ilu rdzeni powinna skorzystać vina				<$cpu>
-n	NAME	Specjalnie nadana nazwa dla próby dokowania.
			Jeśli ten parametr pominięto, nazwa próby dokowania
			zostanie wygenerowana wedle formatki 'E_N_ER_S'.


PRZYKŁADY

Poniższa komenda zadokuje wszystkie ligandy z parametrami
E=40 N=50 ER=50 S=15, korzystając z 6 rdzeni.
Próba dokowania będzie nosić nazwę '40_50_50_15', z formatki.
$ $0 -e 40 -n 50 -r 50 -s 15 -c 6

Poniższa komenda zadokuje wszystkie ligandy z parametrami domyślnymi,
korzystając z 1 rdzenia. Próba dokowania będzie nosić nazwę 'default_params'.
$ $0 -c 1 -n default_params
"
}


#Zczytaj parametry dokowania i nazwę z argumentów do skryptu
while [ "$#" -ge 1 ]; do
	case $1 in
		-h|--help)
			usage
			exit 1
			;;
		-e)
			shift
			e=$1
			;;
		-n)
			shift
			n=$1
			;;
		-r)
			shift
			er=$1
			;;
		-s)
			shift
			size=$1
			;;
		-c)
			shift
			cpu=$1
			;;
		-n)
			shift
			name=$1
			;;
	esac
	shift
done

[ -z "$name" ] && name="${e}_${n}_${er}_${size}"


echo "Nastąpi dokowanie z następującymi parametrami:
Exhaustiveness						$e
Number of modes						$n
Energy ratio						$er
Rozmiar przestrzeni dokowania		$size
Do dokowania wykorzystane zostanie $cpu rdzeni.
Nazwa próby dokowania to $name.
Czy rozpocząć dokowanie? [Y/n]"
read c
[ "$c" = "N" -o "$c" = "n" ] && exit 0


#Kotwica do pętli po folderach
startdir=$PWD

#Pętla po folderach z ligandami
for d in $(find $LIGDIR -mindepth 1 -maxdepth 1 -type d); do

	cd "$d"

	if [ -f "$name/out.pdbqt" ] ; then
		echo "Ligand $d już odbył pomyślną próbę dokowania pod nazwą $name, pomijam."
		cd "$startdir"
		continue
	fi
	if [ ! -f "./$LIGFILE" ] ; then
		echo "Nie znaleziono pliku PDBQT ligandu."
		cd "$startdir"
		continue
	fi

	mkdir -p "$name"

	$VINA \
		--receptor $RECEPTOR --ligand $LIGFILE --cpu $cpu \
		--exhaustiveness $e --num_modes $n --energy_range $er \
		--size_x $size --size_y $size --size_z $size \
		--center_x -16.273 --center_y 31.647 --center_z -26.751 --verbosity 2 \
		--out "$name/out.pdbqt" | tee "$name/out.log"

	cd "$startdir"

done