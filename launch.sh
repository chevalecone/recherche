echo "compiling"
make

echo "Removing previous simulation"
cd output
rm *.vtk
cd ..

echo "Launching simulation"
./lbm.exe


