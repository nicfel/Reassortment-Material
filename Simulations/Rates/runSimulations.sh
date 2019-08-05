files=sim*.xml
for f in $files
do
bsub -J sim -W 4:00 java -jar ~/jar/Coevo.jar $f
done
