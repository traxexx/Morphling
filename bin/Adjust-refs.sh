#!/bin/bash

# set options based on LHMEI directory

while getopts "i:" opt; do
	case $opt in
		i)
			LHMEI_DIR=$OPTARG
			;;
		\?)
			exit 1
			;;
	esac
done

if [ -z $LHMEI_DIR ]; then
	echo -e "\nUsage:"
	echo -e "\t-i\t LHMEI full directory\n"
	exit 1
fi

echo -e "ALU\t$LHMEI_DIR/refs/Alu.fa" > $LHMEI_DIR/refs/MobileElement.list
echo -e "L1\t$LHMEI_DIR/refs/L1.fa" >> $LHMEI_DIR/refs/MobileElement.list
echo -e "SVA\t$LHMEI_DIR/refs/SVA.fa" >> $LHMEI_DIR/refs/MobileElement.list

echo "Finished Creating $LHMEI_DIR/refs/MobileElement.list!"

echo -e "ALU\t$LHMEI_DIR/refs/Alu.bed" > $LHMEI_DIR/refs/MobileElement.coord
echo -e "L1\t$LHMEI_DIR/refs/L1.bed" >> $LHMEI_DIR/refs/MobileElement.coord
echo -e "SVA\t$LHMEI_DIR/refs/SVA.bed" >> $LHMEI_DIR/refs/MobileElement.coord

echo "Finished Creating $LHMEI_DIR/refs/MobileElement.coord!"

$LHMEI_DIR/bin/Lift-over-ctrl-coord.pl -l  $LHMEI_DIR/refs/chr20-slice-MEI.bed -c $LHMEI_DIR/refs/Alu.bed 
$LHMEI_DIR/bin/Lift-over-ctrl-coord.pl -l  $LHMEI_DIR/refs/chr20-slice-MEI.bed -c $LHMEI_DIR/refs/L1.bed 
$LHMEI_DIR/bin/Lift-over-ctrl-coord.pl -l  $LHMEI_DIR/refs/chr20-slice-MEI.bed -c $LHMEI_DIR/refs/SVA.bed 

echo "Finished lift over ctrl bed!"

rm -f $LHMEI_DIR/refs/Alu.bed.unoverlap
rm -f $LHMEI_DIR/refs/L1.bed.unoverlap
rm -f $LHMEI_DIR/refs/SVA.bed.unoverlap

echo -e "ALU\t$LHMEI_DIR/refs/Alu.bed.liftOver" > $LHMEI_DIR/refs/MobileElement.coord.liftOver
echo -e "L1\t$LHMEI_DIR/refs/L1.bed.liftOver" >> $LHMEI_DIR/refs/MobileElement.coord.liftOver
echo -e "SVA\t$LHMEI_DIR/refs/SVA.bed.liftOver" >> $LHMEI_DIR/refs/MobileElement.coord.liftOver

echo "Finished creating $LHMEI_DIR/refs/MobileElement.coord.liftOver!"


