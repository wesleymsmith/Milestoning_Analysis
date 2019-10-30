
trajDir=$1
outDir=$2
trajPrefix=$3

trajSuffix=".traj"
fortBinDir='.'

indexTool='writeindex.exe'
rateTool='writerate.exe'
qtool='qelements.exe'

echo "computing rate matrix and counts for Nij_alpha and Ri_alpha"
escapeMatFile=${outDir}/${trajPrefix}.rate_matrix.dat
rijFile=${outDir}/${trajPrefix}.Rij.dat
nijFile=${outDir}/${trajPrefix}.Nij.dat
tFile=${outDir}/${trajPrefix}.T.dat
rm $escapeMatFile $rijFile $nijFile $tFile 
touch $escapeMatFile $rijFile $nijFile $tFile
for trajFilePath in ${trajDir}/$trajPrefix.*
do 
trajFile=`echo $trajFilePath | sed "s:$trajDir/::g"`
tmpWin=`echo $trajFile | sed "s/$trajPrefix//g" | sed "s/$trajSuffix//g" | sed "s/[.]//g"`
window=`printf "%g" $tmpWin`
echo "working on $trajFile (window $window)"
cp $trajFilePath tmp.traj
${fortBinDir}/$indexTool
cp fort.37 $outDir/${trajPrefix}.${tmpWin}.index_traj.dat
echo $window > fort.5
cat fort.37 >> fort.5
echo $(($window)) > imother.txt
${fortBinDir}/$rateTool
cat fort.38 >> $escapeMatFile
${fortBinDir}/$qtool
cat tmpallN >> $nijFile 
cat tmpallR >> $rijFile
cat tmpallT >> $tFile
done
