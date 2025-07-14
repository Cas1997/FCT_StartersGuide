nEvents=204
inputFolder="\"<Insert your path here>\""
outputFolder="\"<Insert your path here>\""

# root -l -b -q "photonSpectrum.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Run remake manually

# root -l -b -q "photSpectrumPACut.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Run remake & significance manually
# Run signal_over_background manually

root -l -b -q "pointingAngleDistribution.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"

# root -l -b -q "photSpectrumClosestDistance.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Run remake manually

# root -l -b -q "closestDistanceCutAnalysis.cpp(${outputFolder})"

# root -l -b -q "photSpectrumClosestDistanceCut.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Run remake & significance manually
# Run signal_over_background manually

# root -l -b -q "photSpectrumEPIDCut.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Run remake & significance manually
# Run signal_over_background manually

# root -l -b -q "photonOrigin_improved.cpp(${inputFolder}, ${outputFolder}, ${nEvents})" >> "${outputFolder:1:-1}/photonOrigin.txt"
# Manually insert percentages into schematic overview
# Do this for the window too

# root -l -b -q "bremContributions.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Manually run materialbudget plot
# Manually run bremPrediciton.cpp
# Manually run plotPredictionContribution.cpp

# root -l -b -q "chargedParticleSpectrum.cpp(${inputFolder}, ${outputFolder}, ${nEvents})"
# Manually run remake
# Manually insert medium names
