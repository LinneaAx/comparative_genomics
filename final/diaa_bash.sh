
for file in *.fasta
do 
echo $file
python2 -c "import diaa; diaa.compute_diaa('$file')"
done

