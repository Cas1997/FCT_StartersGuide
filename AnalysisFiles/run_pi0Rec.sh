nEvents=204
inputFolder="\"<Insert your path here>\""
outputFolder="\"<Insert your path here>\""

root -l -b -q "pi0_reconstruction.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"

root -l -b -q "pi0_reconstructionPACut.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"

root -l -b -q "pi0_reconstructionClosestDistanceCut.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"

root -l -b -q "pi0_pTRapidityDist.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"