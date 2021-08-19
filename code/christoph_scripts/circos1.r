library(circlize)
circos1 = function(file1, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, countColors = c("#FFFFFFFF", "#0000FFFF")) {
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
	if(length(countColors) < 2) {
		print("Too few countColors specified. Defaulting to blue.")
		countColors = c("#FFFFFFFF", "#0000FFFF")
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
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$Count)
	A$relCount = A$Count / sum

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(A)
		for(index in 1:nrow(A)) {
			cur = cur + A$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		A = A[1:cutoffIndex ,]
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(Count = 0, id = "START", factors = "A", relCount = 0, accumRelCount = 0), A)

	# Combine data sets
	data = dataA

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}

	# Open output file, set parameters and initialize diagramm
	svg(diagrammFileName)
	circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
	circos.initialize(factors = data$factors, x = data$accumRelCount)

	# Add ring with labels, ticks and colored boxes
	circos.track(factors = data$factors, x = data$relCount, y = data$accumRelCount, bg.border = NA,
		panel.fun = function(x, y) {
			circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(4), fileAliases[CELL_META$sector.numeric.index], niceFacing = TRUE)
			circos.lines(c(CELL_META$xlim[1], CELL_META$xlim[2]), c((CELL_META$ylim[1] + CELL_META$ycenter) * 0.5, (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
			for(index in 1:length(x)) {
				circos.lines(c(y[index], y[index]), c(CELL_META$ylim[1], (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
			}
			n = length(x) - 1
			circos.rect(y[2:(n + 1)] - x[2:(n + 1)], rep(CELL_META$ycenter, n), y[2:(n + 1)], rep(CELL_META$ylim[2], n), col = countRamp(x[2:(n + 1)]), border = NA)
		})

	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos1.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			cufoff = c(cutoff),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount))
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