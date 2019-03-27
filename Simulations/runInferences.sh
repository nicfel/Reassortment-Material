files=inf*.xml
for f in $files
do
bsub -J sim -W 24:00 -R "rusage[mem=8000]" java -jar ~/jar/Coevo.jar $f
done
