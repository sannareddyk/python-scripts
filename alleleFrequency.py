import sys
import random
import os
import re

###############################################################
### usage Section

if len(sys.argv) != 4:
	print("usage: python parseVCF.py vcfFileName nSamples outfileName")
	sys.exit()

###############################################################
#### Inputs


vcfFileName = sys.argv[1]
nSamples = int(sys.argv[2])
outfileName = sys.argv[3]

#### Open input File and start processing
vcfFile = open(vcfFileName)
outfile = open(outfileName,'w')
fieldDict = {}

#### Extra Variables
#sumHeterozygosity = 0
#numVariants = 0

##### For each line in vcf file, 
# 1. If the line is a commented line ( '##' ) , no processing needed
# 2. If the line is a header ( '#' ) , create a field dictionary
# 3. If not 1,2 , then parse the line and process data
# 4. 
for line in vcfFile:

	if line.startswith('##'):
		## This is just a comment line. So just write it to output
		outfile.write(line)
	elif line.startswith('#'):
		## This is a header line. Get the field dict from here.
		outfile.write(line)
		line = line.strip()
		line = line.lstrip('#')
		line = line.split("\t")
		counter = 0
		for f in line:
			fieldDict[f] = counter
			counter += 1
	else:

		### Actual data line, process it here.
		formatColumn = 0
		referenceColumnNum = 0
		alternateColumnNum = 0

		if 'FORMAT' in fieldDict:
			formatColumn = fieldDict['FORMAT']
		if 'REF' in fieldDict:
			referenceColumnNum = fieldDict['REF']
		if 'ALT' in fieldDict:
			alternateColumnNum = fieldDict['ALT']

		if (formatColumn != 0) and (referenceColumnNum != 0) and (alternateColumnNum != 0):

			## start calculations here.
			resultLine=line.strip()
			line = line.strip().split('\t')

			#format = line[formatColumn]
			refAllele = line[referenceColumnNum]
			altAllele = line[alternateColumnNum]

                        if altAllele in ['A','T','G','C']:

				##### If the altAllele is A or T or G or C, do the following.
				tempDict = {}
				sampleNum = 1

				## Read the data into tempDict
				for item in line[formatColumn+1:formatColumn+nSamples+1]:
					tempDict[sampleNum] = item
					sampleNum += 1

				#### Define the genotype dict and calculate heterozygosity
				genotypeDict={}
				for sampleNum in tempDict:

					datum=tempDict[sampleNum]
					datum=datum.strip().split(':')

					GT = datum[0]
					PL = datum[1]
					DP = datum[2]
					SP = datum[3]
					GQ = datum[4]

					## assign genotype
					genotype = '--'
					if GT == '0/0':
						genotype='homRef'
					elif GT == '0/1':
						genotype='het'
					elif GT =='1/1':
						genotype='homAlt'

					genotypeDict[sampleNum] = []

					if genotype=='homRef':
						genotypeDict[sampleNum].append(refAllele)
					elif genotype=='homAlt':
						genotypeDict[sampleNum].append(altAllele)
					elif genotype=='het':            
						tempArray = [refAllele,altAllele]
						randomAllele = random.sample(tempArray,1)[0]
						genotypeDict[sampleNum].append(randomAllele)

					PL = PL.strip().split(',')
					homRef_PL = int(PL[0])
					het_PL = int(PL[1])
					homAlt_PL = int(PL[2])

					isBadCall = 0
					## criteria for bad call
					if ( homRef_PL + het_PL + homAlt_PL  ) == 0 :
						isBadCall = 1
					if int(DP) == 0: 
						isBadCall = 1
					if int(GQ) < 20:
						isBadCall = 1

					genotypeDict[sampleNum].append(isBadCall)

				## Compute heterozygosity for the first 18 samples
				n = 0
				i = 0

				for sampleNum in genotypeDict:
					if sampleNum<=18:
						isBadCall = genotypeDict[sampleNum][1]
						genotype = genotypeDict[sampleNum][0]

						if isBadCall == 0:
							n += 1
							if genotype == altAllele:
								i += 1

						resultLine += '\t' + genotype
						resultLine += '\t' + str(isBadCall)

				if n != 0:
					AF1 = i * 1.0 / n
					AF1 = round(AF1,4)
					heterozygosity =  2.0 * i * 1.0/n * (n-i) * 1.0/n
					resultLine += '\t' + str(i) + '\t' + str(n) + '\t' + str(AF1) + "\t" + heterozygosity
				elif n == 0:
					resultLine += '\t' + str(i) + '\t' + str(n) + '\t' + "--" + "\t" + '--'

				## Compute heterozygosity for the next 18 samples
				n=0
				i=0

				for sampleNum in genotypeDict:
					if sampleNum>18 and sampleNum<=36:
						isBadCall = genotypeDict[sampleNum][1]
						genotype = genotypeDict[sampleNum][0]

						if isBadCall == 0:
							n += 1
							if genotype == altAllele:
								i += 1

						resultLine += '\t' + genotype
						resultLine += '\t' + str(isBadCall)

				if n!=0:
					AF2 = i * 1.0 / n
					AF2 = round(AF2,4)
					heterozygosity =  2.0 * i * 1.0/n * (n-i) * 1.0/n     
					resultLine += '\t' + str(i) + '\t' + str(n) + '\t' + str(AF2) + '\t' + str(heterozygosity)
				elif n == 0:
					resultLine += '\t' + str(i) + '\t' + str(n) + '\t' + "--"+ '\t' + "--"

			## Write the result line to output file. 
			outfile.write(resultLine + '\n')

			## If the altAllele is not in A, T, G, or C, do the following
                        elif altAllele == ".":

			    ## Read data to temp dict    
			    tempDict = {}
			    sampleNum = 1
			    for item in line[formatColumn+1:formatColumn+nSamples+1]:
			        tempDict[sampleNum] = item
				sampleNum += 1

			    ## create genotype dict from temp Dict
			    genotypeDict={}
			    i=0
			    n=0
			    AF1 = 0
			    AF2 = 0             
			    heterozygosity = 0

			    ## for samples 1 to 18 do the following
			    for sampleNum in tempDict:
			        if sampleNum<=18:
				    datum=tempDict[sampleNum]
				    datum=datum.strip().split(':')

				    PL = datum[0]
				    DP = datum[1]
				    SP = datum[2]

				    ## assign genotype
				    genotype = '--'                        
				    genotypeDict[sampleNum] = []
				    genotypeDict[sampleNum].append(refAllele)

				    isBadCall=1
				    genotypeDict[sampleNum].append(isBadCall)

				    isBadCall = genotypeDict[sampleNum][1]
				    genotype = genotypeDict[sampleNum][0]

				    resultLine += '\t' + genotype
				    resultLine += '\t' + str(isBadCall)

				## Append the summary of the first 18 samples        
				resultLine += '\t' + str(i) + '\t' + str(n) + '\t' + str(AF1) + '\t' + str(heterozygosity)

				## For samples 19 to 86 , do the following
				for sampleNum in tempDict:
				    if sampleNum>18 and sampleNum<=36:
                                        datum=tempDict[sampleNum]
					datum=datum.strip().split(':')

					PL = datum[0]
					DP = datum[1]
					SP = datum[2]

					## assign genotype
					genotype = '--'

					genotypeDict[sampleNum] = []
					genotypeDict[sampleNum].append(refAllele)

					isBadCall=1
					genotypeDict[sampleNum].append(isBadCall)
						
					isBadCall = genotypeDict[sampleNum][1]
					genotype = genotypeDict[sampleNum][0]

					resultLine += '\t' + genotype
					resultLine += '\t' + str(isBadCall)

				### write the summary of 19 to 36 samples            
				resultLine += '\t' + str(i) + '\t' + str(n) + '\t'+ str(AF2) + '\t' + str(heterozygosity)

				## Write the final output
				outfile.write(resultLine+ '\n')
                    
#print('sumHeterozygosity: '+str(sumHeterozygosity))
#print('numVariants: '+str(numVariants))
outfile.close()
vcfFile.close()
