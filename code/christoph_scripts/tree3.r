library(ggplot2)
library(treemapify)
tree3 = function(file1, file2, file3, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, colors = c("#FF0000", "#00FF00", "#0000FF", "#5A5A5A")) {
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
	if(length(colors) < 4) {
		print("Too few colors specified. Defaulting to red, green, blue and grey.")
		colors = c("#FF0000", "#00FF00", "#0000FF", "#5A5A5A")
	}

	# Helper functions
	readData = function(file) {
		data = read.table(file, header = TRUE, sep = "\t", colClasses = c("integer", "character", "character", "character", "character", "numeric"))
		data = cbind(data, data.frame(id = paste(data$CDR3.nucleotide.sequence, data$V.segments, data$J.segments, sep="")))
		subset(data, select = -c(CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, V.segments, J.segments, Percentage))
	}
	getStatus = function(data) {
		temp = data.frame(one = unlist(data[4:6]), two = c("A", "B", "C"))
		statusLetters = paste(temp[temp$one > 0, "two"], collapse = "")
		paste(nchar(statusLetters), statusLetters, sep = "")
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
	diagrammFileName = paste(saveFolder, diagrammName, "_tree.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = paste("A:", fileAliases[1])
	A$linkIndexA = 1
	A$linkIndexB = 0
	A$linkIndexC = 0

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
	B$linkIndexA = 0
	B$linkIndexB = 1
	B$linkIndexC = 0

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
	C$linkIndexA = 0
	C$linkIndexB = 0
	C$linkIndexC = 1

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

	# Compute links between A and B
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
			}
		}
	}

	# Compute links between A and C
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(C)) {
			if(A$id[index_1] == C$id[index_2]){
				A$linkIndexC[index_1] = index_2
				C$linkIndexA[index_2] = index_1
			}
		}
	}

	# Compute links between B and C
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(C)) {
			if(B$id[index_1] == C$id[index_2]){
				B$linkIndexC[index_1] = index_2
				C$linkIndexB[index_2] = index_1
			}
		}
	}

	# Set status
	for(index in 1:nrow(A)) {
		A$status[index] = getStatus(A[index,])
	}

	# Set status
	for(index in 1:nrow(B)) {
		B$status[index] = getStatus(B[index,])
	}

	# Set status
	for(index in 1:nrow(C)) {
		C$status[index] = getStatus(C[index,])
	}

	# Combine data sets
	data = rbind(A, B, C)

	# Remove colors corresponding to not existing share statuses
	status = c("1A", "1B", "1C", "2AB", "2AC", "2BC", "3ABC")
	colors = append(c("#FFFFFF", "#FFFFFF", "#FFFFFF"), colors)
	temp = data.frame(exists = logical(), color = character())
	for(index in 1:length(status)) {
		temp = rbind(temp, data.frame(exists = c(nrow(data[data$status == status[index],]) > 0), color = c(colors[index])))
	}
	colors = as.vector(unlist(temp[temp$exists == TRUE, "color"])) # In R studio as.vector(unlist()) is not necessary but in the command line it is

	# Create plot
	p <- ggplot(data, aes(area = Count, fill = status)) + 
		geom_treemap(colour = "#000000") + 
		scale_fill_manual(values = colors) +
		theme(aspect.ratio = 1, strip.background = element_rect(fill="#FFFFFF")) +
		facet_wrap( ~ factors)

	# Save plot
	svg(diagrammFileName, width = 12, height = 4)
	print(p)
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "tree3.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			cufoff = c(cutoff),
			colors1 = c(colors[1]),
			colors2 = c(colors[2]),
			colors3 = c(colors[3]),
			colors4 = c(colors[4]))
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