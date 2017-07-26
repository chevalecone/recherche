echo "compiling"
make

echo "Deleting previous results"
rm output/*

echo "Launching simulation anew"
./lbm.exe

