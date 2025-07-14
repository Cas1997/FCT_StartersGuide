# To run cherenkov detector simulations
sipmFileBase="\"/input_SiPM_pdes/"

outputFolderHe="\"\""
paramFileBaseHe="\"/paramFiles/helium/"
outputTxtBaseHe="/helium/outputTxt"

declare -a paramArr=("hamamatsu_3050cs" "hamamatsu_3075cs")
for i in "${paramArr[@]}"
do
	for ((param=1; param<7; param++));
	do
		paramFileName=${paramFileBaseHe}"$i""_param""$param"".cfg\""
		sipmFileName=${sipmFileBase}"$i""_pde.cfg\""
		outputTxtName=${outputTxtBaseHe}"_""$i""_param""$param"".txt"
		fileNameAddition="$i""_param""$param"
		root -l -b -q "CherenkovDetector/cdEventVeto.cpp(${inputFolder}, ${outputFolderHe}, ${paramFileName}, ${sipmFileName}, \"${fileNameAddition}\", ${nEvents})" >> $outputTxtName
	done
done

for ((param=1; param<7; param++));
do
	paramFileName=${paramFileBaseHe}"FBK_param""$param"".cfg\""
	sipmFileName=${sipmFileBase}"FBK_pde_3V.cfg\""
	outputTxtName=${outputTxtBaseHe}"_FBK_param""$param"".txt"
	fileNameAddition="FBK_param""$param"
	root -l -b -q "CherenkovDetector/cdEventVeto.cpp(${inputFolder}, ${outputFolderHe}, ${paramFileName}, ${sipmFileName}, \"${fileNameAddition}\", ${nEvents})" >> $outputTxtName
done

outputFolderNe="\"/neon\""
paramFileBaseNe="\"/paramFiles/neon/"
outputTxtBaseNe="/neon/outputTxt"

for i in "${paramArr[@]}"
do
	for ((param=1; param<7; param++));
	do
		paramFileName=${paramFileBaseNe}"$i""_param""$param"".cfg\""
		sipmFileName=${sipmFileBase}"$i""_pde.cfg\""
		outputTxtName=${outputTxtBaseNe}"_""$i""_param""$param"".txt"
		fileNameAddition="$i""_param""$param"
		root -l -b -q "CherenkovDetector/cdEventVeto.cpp(${inputFolder}, ${outputFolderNe}, ${paramFileName}, ${sipmFileName}, \"${fileNameAddition}\", ${nEvents})" >> $outputTxtName
	done
done

for ((param=1; param<7; param++));
do
	paramFileName=${paramFileBaseNe}"FBK_param""$param"".cfg\""
	sipmFileName=${sipmFileBase}"FBK_pde_3V.cfg\""
	outputTxtName=${outputTxtBaseNe}"_FBK_param""$param"".txt"
	fileNameAddition="FBK_param""$param"
	root -l -b -q "CherenkovDetector/cdEventVeto.cpp(${inputFolder}, ${outputFolderNe}, ${paramFileName}, ${sipmFileName}, \"${fileNameAddition}\", ${nEvents})" >> $outputTxtName
done

# Manually run remake
