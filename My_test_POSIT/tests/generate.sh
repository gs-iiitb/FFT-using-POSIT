
if [ "$#" != "1" ]
then
    echo -e "Please give No.of Times to run.\n"
    exit
fi

rm smoke_tests*.txt
 
echo -e "Generating $1 Test cases...\n"
END=$1

for i in $(seq 1 $END)
do 
    ./smoke_randoms > smoke_tests_${i}.txt
done

echo -e "DONE: Generating.\n"

echo -e "Now Separating Operations...\n"


for i in $(seq 1 $END)
do
   echo -e "  Extracting Div Operation from: smoke_test_$i\n"
   sed -n -e '306,404p' smoke_tests_${i}.txt >> smoke_tests_rec.txt
done
echo -e "DONE: DIV.\n"


