library(ggplot2)
library(treemapify)
track3 = function(file1, file2, file3, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, numToTrack = 5, baseFileId = 1, baseFile = NULL, colors = rainbow(5)) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 3) {
		print("To few file aliases specified. Defaulting to file names.")
		fileAliases = NULL
	}
	if(cutoff < 0.0) {
		pirnt("Cutoff value is below zero. Defaulting to 1.0.")
		cutoff = 1.0
	}
	if(cutoff > 1.0) {
		print("Cutoff value is above 1.0. Defaulting to 1.0.")
		cutoff = 1.0
	}
	if(numToTrack < 1) {
		print("Number of clones track is bolow zero. Defaulting to 5.")
		numToTrack = 5
	}
	if(baseFileId < 0) {
		print("Base file id is below zero. Defaulting to 1.")
		baseFileId = 1
	}
	if(baseFileId == 0 && is.null(baseFile)) {
		print("No baseFile is specified. Defaulting to file1.")
		baseFileId = 1
	}
	if(baseFileId > 3) {
		print("Base file id is above 3. Defaulting to 1.")
		baseFileId = 1
	}

	# Helper functions
	readData = function(file) {
		data = read.table(file, header = TRUE, sep = "\t", colClasses = c("integer", "character", "character", "character", "character", "numeric"))
		data = cbind(data, data.frame(id = paste(data$CDR3.nucleotide.sequence, data$V.segments, data$J.segments, sep="")))
		subset(data, select = -c(CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, V.segments, J.segments, Percentage))
	}

	# Extract file names
	file1Parts = unlist(strsplit(file1, "/"))
	file1Name = unlist(strsplit(file1Parts[length(file1Parts)], "[.]"))[1]
	file2Parts = unlist(strsplit(file2, "/"))
	file2Name = unlist(strsplit(file2Parts[length(file2Parts)], "[.]"))[1]
	file3Parts = unlist(strsplit(file3, "/"))
	file3Name = unlist(strsplit(file3Parts[length(file3Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name, file2Name, file3Name)
	}

	# Create save folder name if not provided
	if(is.null(saveFolder)) {
		if(length(file1Parts) > 1) {
			saveFolder = paste(file1Parts[1:length(file1Parts) - 1], collapse = "/")
			saveFolder = paste(saveFolder, "/", sep = "")
		} else {
			saveFolder = ""
		}
	}

	# Create diagramm name
	diagrammName = paste(file1Name, file2Name, file3Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_track.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = paste("A:", fileAliases[1])
	A$status = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(A)
		sum = sum(A$Count)
		for(index in 1:nrow(A)) {
			cur = cur + (A$Count[index] / sum)
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		A = A[1:cutoffIndex ,]
	}

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = paste("B:", fileAliases[2])
	B$status = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		sum = sum(B$Count)
		for(index in 1:nrow(B)) {
			cur = cur + (B$Count[index] / sum)
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Read data from file3 and add additional columns
	C = readData(file3)
	C$factors = paste("C:", fileAliases[3])
	C$status = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(C)
		sum = sum(C$Count)
		for(index in 1:nrow(C)) {
			cur = cur + (C$Count[index] / sum)
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		C = C[1:cutoffIndex ,]
	}

	# Extract clones to track
	if(baseFileId == 0) {
		baseFileParts = unlist(strsplit(baseFile, "/"))
		baseFileName = unlist(strsplit(baseFileParts[length(baseFileParts)], "[.]"))[1]
		baseFileData = readData(baseFile)
		numToTrack = nrow(baseFileData)
		toTrack = baseFileData[1:numToTrack,"id"]
	}
	if(baseFileId == 1) {
		baseFileName = NA
		toTrack = A[1:numToTrack,"id"]
	}
	if(baseFileId == 2) {
		baseFileName = NA
		toTrack = B[1:numToTrack,"id"]
	}
	if(baseFileId == 3) {
		baseFileName = NA
		toTrack = C[1:numToTrack,"id"]
	}

	# Update colors
	if(length(colors) != numToTrack) {
		print("Number of supplied colors does not match number of tracked clones. Defaulting to rainbow colors.")
		colors = rainbow(numToTrack)
	}

	# Track clones in A
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:length(toTrack)) {
			if(A$id[index_1] == toTrack[index_2]){
				A$status[index_1] = index_2
			}
		}
	}

	# Track clones in B
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:length(toTrack)) {
			if(B$id[index_1] == toTrack[index_2]){
				B$status[index_1] = index_2
			}
		}
	}

	# Track clones in C
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:length(toTrack)) {
			if(C$id[index_1] == toTrack[index_2]){
				C$status[index_1] = index_2
			}
		}
	}

	# Combine data sets
	data = rbind(A, B, C)

	# Set status
	data$status = formatC(data$status, width = floor(log10(numToTrack)) + 1, format = "d", flag = "0")

	# Add color for non tracked clones
	colors = append(c("#FFFFFF"), colors)

	# Helper function that lets us hide the legend
	legend <- function(flag) {
		if (flag) {
			theme(legend.position = "none")
		}
	}

	# Create plot
	p <- ggplot(data, aes(area = Count, fill = status)) + 
		geom_treemap(colour = "#000000") + 
		scale_fill_manual(values = colors, name = "Tracked clones") +
		theme(aspect.ratio = 1, strip.background = element_rect(fill="#FFFFFF")) +
		legend(numToTrack > 50) +
		facet_wrap( ~ factors)

	# Save plot
	svg(diagrammFileName, width = 12, height = 4)
	print(p)
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "track3.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			cufoff = c(cutoff),
			numToTrack = c(numToTrack),
			baseFileId = c(baseFileId),
			baseFileName = c(baseFileName))
	if(file.exists(parameterFileName) == TRUE) {
		parameters = read.csv(parameterFileName)
		if(nrow(parameters[parameters$diagrammName == diagrammName, ]) == 0) {
			parameters = rbind(parameters, parameter)
		} else {
			previousIndex = match(c(diagrammName), parameters$diagrammName)
			parameters[previousIndex,] = parameter
		}
		write.csv(parameters, parameterFileName, row.names = FALSE)
	} else {
		write.csv(parameter, parameterFileName, row.names = FALSE)
	}
}
