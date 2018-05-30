#cd ../splitoutput/last/ 

#genome files
#cd ./genomes/
#for file in *.fa.txt
#do 
#echo $file
#python2 -c 'import grand_finale.py; ORF_finder(\"$file\")'
#python2 -c 'import grand_finale.py; compute_gc(\"$file\")'
#python2 -c 'import grand_finale.py; compute_dinucleo(\"$file\")'
#python2 -c 'import grand_finale.py; compute_dinucleo(\"$file\")'
#done



#cd ./genomes/
#for file in *.fa.txt
#do 
#    ls | tr "\n" " " | $file 
#done 

#cd ./genomes/
python2 -c "import grand_finale; grand_finale.distance_matrix('$file')"

#python2 grand_finale.py 03.fa.txt 28.fa.txt 43.fa.txt 48.fa.txt 50.fa.txt

for file in *.fa.txt
do 
echo $file 
python2 grand_finale.py $file
done

for file in *.fasta
do 
echo $file
python2 -c "import grand_finale; grand_finale.compute_aa('$file')"
done
