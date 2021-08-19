library(ggplot2)
library(treemapify)
tree1 = function(file1, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 1) {
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

	# Helper functions
	readData = function(file) {
		data = read.table(file, header = TRUE, sep = "\t", colClasses = c("integer", "character", "character", "character", "character", "numeric"))
		data = cbind(data, data.frame(id = paste(data$CDR3.nucleotide.sequence, data$V.segments, data$J.segments, sep="")))
		subset(data, select = -c(CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, V.segments, J.segments, Percentage))
	}

	# Extract file names
	file1Parts = unlist(strsplit(file1, "/"))
	file1Name = unlist(strsplit(file1Parts[length(file1Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name)
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
	diagrammName = paste(file1Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_tree.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = paste("A:", fileAliases[1])
	A$status = "1A"

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

	# Combine data sets
	data = A

	# Create plot
	p <- ggplot(data, aes(area = Count, fill = status)) + 
		geom_treemap(colour = "#000000") + 
		scale_fill_manual(values = c("#FFFFFF")) +
		theme(aspect.ratio = 1, strip.background = element_rect(fill="#FFFFFF")) +
		facet_wrap( ~ factors)

	# Save plot
	svg(diagrammFileName, width = 4, height = 4)
	print(p)
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "tree1.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			cufoff = c(cutoff))
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