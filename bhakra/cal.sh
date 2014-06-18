rm ~/projects/pakistan/results.txt
touch ~/projects/pakistan/results.txt
#cp ~/projects/pakistan/saved2.gas ~/projects/pakistan/saved3.gas
#cp ~/projects/pakistan/saved.gas ~/projects/pakistan/saved2.gas
cp saved2.gas saved3.gas
cp saved.gas saved2.gas
for n in {1..10}
do
    echo $n
    nice ./model >> ~/projects/pakistan/results.txt
done
