library(circlize)

# Helper functions
readData = function(file) {
    data = read.csv(file, header = TRUE, sep = "\t")
    data = cbind(data, data.frame(id = paste(data$nSeqCDR3, data$bestVHit, data$bestJHit, sep="_")))
    subset(data, select = c(id, cloneCount))
}

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
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum

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
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, accumRelCount = 0), A)

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

circos2 = function(file1, file2, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080")) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 2) {
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
	if(length(linkColors) < 1) {
		print("Too few linkColors specified. Defaulting to red.")
		linkColors = c("#FF000080", "#FF000080")
	}
	if(length(linkColors) == 1) {
		linkColors = c(linkColors[1], linkColors[1])
	}

	getStatus = function(data) {
		temp = data.frame(one = unlist(data[5:6]), two = c("A", "B"))
		statusLetters = paste(temp[temp$one > 0, "two"], collapse = "")
		paste(nchar(statusLetters), statusLetters, sep = "")
	}

	# Extract file names
	file1Parts = unlist(strsplit(file1, "/"))
	file1Name = unlist(strsplit(file1Parts[length(file1Parts)], "[.]"))[1]
	file2Parts = unlist(strsplit(file2, "/"))
	file2Name = unlist(strsplit(file2Parts[length(file2Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name, file2Name)
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
	diagrammName = paste(file1Name, file2Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum
	A$linkIndexA = 1
	A$linkIndexB = 0

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

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = "B"
	sum = sum(B$cloneCount)
	B$relCount = B$cloneCount / sum
	B$linkIndexA = 0
	B$linkIndexB = 1

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		for(index in 1:nrow(B)) {
			cur = cur + B$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Compute links between A and B
	maxLinkSize = 0
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + B$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
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

	if(sort == TRUE) {
		# Sort by status
		A = A[order(A$status),]

		# Update possible links
		for(index in 1:nrow(A)) {
			if(A$linkIndexB[index] > 0) {
				B$linkIndexA[A$linkIndexB[index]] = index
			}
		}

		# Sort by status
		B = B[order(B$status),]

		# Update possible links
		for(index in 1:nrow(B)) {
			if(B$linkIndexA[index] > 0) {
				A$linkIndexB[B$linkIndexA[index]] = index
			}
		}
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(B)) {
		cur = cur + B$relCount[index]
		B$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, status = "1A", accumRelCount = 0), A)
	dataB = rbind(data.frame(cloneCount = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, status = "1B", accumRelCount = 0), B)

	# Combine data sets
	data = rbind(dataA, dataB)

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	
	# Define link color ramp for each type of link
	linkRampAB = function(x) {
		ramp = colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
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

	# Add links
	for(index in 1:nrow(A)) {
		if(A$linkIndexB[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRampAB(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
		}
	}

	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos2.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			cufoff = c(cutoff),
			sort = c(sort),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount),
			linkColorsAB1 = c(linkColors[1]),
			linkColorsAB2 = c(linkColors[2]),
			maxLinkSize = c(maxLinkSize))
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

circos3 = function(file1, file2, file3, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080", "#FF000080", "#FF000080"), showLinks = c(TRUE, TRUE, TRUE)) {
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
	if(length(countColors) < 2) {
		print("Too few countColors specified. Defaulting to blue.")
		countColors = c("#FFFFFFFF", "#0000FFFF")
	}
	if(length(linkColors) < 3) {
		print("Too few linkColors specified. Defaulting to red.")
		linkColors = rep(c("#FF000080", "#FF000080", "#FF000080"), each = 2)
	}
	if(length(linkColors) == 3) {
		linkColors = rep(linkColors, each = 2)
	}
	if(length(linkColors) > 3 && length(linkColors) < 6) {
		print("Too few linkColors for link color ramp specified. Defaulting to red.")
		linkColors = rep(c("#FF000080", "#FF000080", "#FF000080"), each = 2)
	}
	if(length(showLinks) < 3) {
		print("Too few values in showLinks specified. Defaulting to showing all links.")
		showLinks = c(TRUE, TRUE, TRUE)
	}

	getStatus = function(data) {
		temp = data.frame(one = unlist(data[5:7]), two = c("A", "B", "C"))
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
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum
	A$linkIndexA = 1
	A$linkIndexB = 0
	A$linkIndexC = 0

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

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = "B"
	sum = sum(B$cloneCount)
	B$relCount = B$cloneCount / sum
	B$linkIndexA = 0
	B$linkIndexB = 1
	B$linkIndexC = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		for(index in 1:nrow(B)) {
			cur = cur + B$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Read data from file3 and add additional columns
	C = readData(file3)
	C$factors = "C"
	sum = sum(C$cloneCount)
	C$relCount = C$cloneCount / sum
	C$linkIndexA = 0
	C$linkIndexB = 0
	C$linkIndexC = 1

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(C)
		for(index in 1:nrow(C)) {
			cur = cur + C$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		C = C[1:cutoffIndex ,]
	}

	# Compute links between A and B
	maxLinkSize = 0
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + B$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and C
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(C)) {
			if(A$id[index_1] == C$id[index_2]){
				A$linkIndexC[index_1] = index_2
				C$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and C
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(C)) {
			if(B$id[index_1] == C$id[index_2]){
				B$linkIndexC[index_1] = index_2
				C$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
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

	if(sort == TRUE) {
		# Sort by status
		A = A[order(A$status),]

		# Update possible links
		for(index in 1:nrow(A)) {
			if(A$linkIndexB[index] > 0) {
				B$linkIndexA[A$linkIndexB[index]] = index
			}
			if(A$linkIndexC[index] > 0) {
				C$linkIndexA[A$linkIndexC[index]] = index
			}
		}

		# Sort by status
		B = B[order(B$status),]

		# Update possible links
		for(index in 1:nrow(B)) {
			if(B$linkIndexA[index] > 0) {
				A$linkIndexB[B$linkIndexA[index]] = index
			}
			if(B$linkIndexC[index] > 0) {
				C$linkIndexB[B$linkIndexC[index]] = index
			}
		}

		# Sort by status
		C = C[order(C$status),]

		# Update possible links
		for(index in 1:nrow(C)) {
			if(C$linkIndexA[index] > 0) {
				A$linkIndexC[C$linkIndexA[index]] = index
			}
			if(C$linkIndexB[index] > 0) {
				B$linkIndexC[C$linkIndexB[index]] = index
			}
		}
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(B)) {
		cur = cur + B$relCount[index]
		B$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(C)) {
		cur = cur + C$relCount[index]
		C$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1A", accumRelCount = 0), A)
	dataB = rbind(data.frame(cloneCount = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1B", accumRelCount = 0), B)
	dataC = rbind(data.frame(cloneCount = 0, id = "START", factors = "C", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1C", accumRelCount = 0), C)

	# Combine data sets
	data = rbind(dataA, dataB, dataC)

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}

	# Define link color ramp for each type of link
	linkRampAB = function(x) {
		ramp = colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampAC = function(x) {
		ramp = colorRamp(c(linkColors[3], linkColors[4]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampBC = function(x) {
		ramp = colorRamp(c(linkColors[5], linkColors[6]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
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

	# Add links
	for(index in 1:nrow(A)) {
		if(A$linkIndexB[index] > 0 && showLinks[1] == TRUE) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRampAB(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
		}
		if(A$linkIndexC[index] > 0 && showLinks[2] == TRUE) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "C", c(C$accumRelCount[A$linkIndexC[index]] - C$relCount[A$linkIndexC[index]], C$accumRelCount[A$linkIndexC[index]]), col = linkRampAC(A$relCount[index] + C$relCount[A$linkIndexC[index]]), border = NA)
		}
	}
	for(index in 1:nrow(B)) {
		if(B$linkIndexC[index] > 0 && showLinks[3] == TRUE) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "C", c(C$accumRelCount[B$linkIndexC[index]] - C$relCount[B$linkIndexC[index]], C$accumRelCount[B$linkIndexC[index]]), col = linkRampBC(B$relCount[index] + C$relCount[B$linkIndexC[index]]), border = NA)
		}
	}

	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos3.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			cufoff = c(cutoff),
			sort = c(sort),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount),
			linkColorsAB1 = c(linkColors[1]),
			linkColorsAB2 = c(linkColors[2]),
			linkColorsAC1 = c(linkColors[3]),
			linkColorsAC2 = c(linkColors[4]),
			linkColorsBC1 = c(linkColors[5]),
			linkColorsBC2 = c(linkColors[6]),
			maxLinkSize = c(maxLinkSize),
			showLinksAB = c(showLinks[1]),
			showLinksAC = c(showLinks[2]),
			showLinksBC = c(showLinks[3]))
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

circos4 = function(file1, file2, file3, file4, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080"), showLinks = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 4) {
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
	if(length(linkColors) < 6) {
		print("Too few linkColors specified. Defaulting to red.")
		linkColors = rep(c("#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080"), each = 2)
	}
	if(length(linkColors) == 6) {
		linkColors = rep(linkColors, each = 2)
	}
	if(length(linkColors) > 6 && length(linkColors) < 12) {
		print("Too few linkColors for link color ramp specified. Defaulting to red.")
		linkColors = rep(c("#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080", "#FF000080"), each = 2)
	}
	if(length(showLinks) < 6) {
		print("Too few values in showLinks specified. Defaulting to showing all links.")
		showLinks = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
	}

	getStatus = function(data) {
		temp = data.frame(one = unlist(data[5:8]), two = c("A", "B", "C", "D"))
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
	file4Parts = unlist(strsplit(file4, "/"))
	file4Name = unlist(strsplit(file4Parts[length(file4Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name, file2Name, file3Name, file4Name)
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
	diagrammName = paste(file1Name, file2Name, file3Name, file4Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum
	A$linkIndexA = 1
	A$linkIndexB = 0
	A$linkIndexC = 0
	A$linkIndexD = 0

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

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = "B"
	sum = sum(B$cloneCount)
	B$relCount = B$cloneCount / sum
	B$linkIndexA = 0
	B$linkIndexB = 1
	B$linkIndexC = 0
	B$linkIndexD = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		for(index in 1:nrow(B)) {
			cur = cur + B$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Read data from file3 and add additional columns
	C = readData(file3)
	C$factors = "C"
	sum = sum(C$cloneCount)
	C$relCount = C$cloneCount / sum
	C$linkIndexA = 0
	C$linkIndexB = 0
	C$linkIndexC = 1
	C$linkIndexD = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(C)
		for(index in 1:nrow(C)) {
			cur = cur + C$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		C = C[1:cutoffIndex ,]
	}

	# Read data from file4 and add additional columns
	D = readData(file4)
	D$factors = "D"
	sum = sum(D$cloneCount)
	D$relCount = D$cloneCount / sum
	D$linkIndexA = 0
	D$linkIndexB = 0
	D$linkIndexC = 0
	D$linkIndexD = 1

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(D)
		for(index in 1:nrow(D)) {
			cur = cur + D$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		D = D[1:cutoffIndex ,]
	}

	# Compute links between A and B
	maxLinkSize = 0
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + B$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between file1 and C
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(C)) {
			if(A$id[index_1] == C$id[index_2]){
				A$linkIndexC[index_1] = index_2
				C$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and D
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(D)) {
			if(A$id[index_1] == D$id[index_2]){
				A$linkIndexD[index_1] = index_2
				D$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and C
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(C)) {
			if(B$id[index_1] == C$id[index_2]){
				B$linkIndexC[index_1] = index_2
				C$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and D
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(D)) {
			if(B$id[index_1] == D$id[index_2]){
				B$linkIndexD[index_1] = index_2
				D$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and D
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(D)) {
			if(C$id[index_1] == D$id[index_2]){
				C$linkIndexD[index_1] = index_2
				D$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
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

	# Set status
	for(index in 1:nrow(D)) {
		D$status[index] = getStatus(D[index,])
	}

	if(sort == TRUE) {
		# Sort by status
		A = A[order(A$status),]

		# Update possible links
		for(index in 1:nrow(A)) {
			if(A$linkIndexB[index] > 0) {
				B$linkIndexA[A$linkIndexB[index]] = index
			}
			if(A$linkIndexC[index] > 0) {
				C$linkIndexA[A$linkIndexC[index]] = index
			}
			if(A$linkIndexD[index] > 0) {
				D$linkIndexA[A$linkIndexD[index]] = index
			}
		}

		# Sort by status
		B = B[order(B$status),]

		# Update possible links
		for(index in 1:nrow(B)) {
			if(B$linkIndexA[index] > 0) {
				A$linkIndexB[B$linkIndexA[index]] = index
			}
			if(B$linkIndexC[index] > 0) {
				C$linkIndexB[B$linkIndexC[index]] = index
			}
			if(B$linkIndexD[index] > 0) {
				D$linkIndexB[B$linkIndexD[index]] = index
			}
		}

		# Sort by status
		C = C[order(C$status),]

		# Update possible links
		for(index in 1:nrow(C)) {
			if(C$linkIndexA[index] > 0) {
				A$linkIndexC[C$linkIndexA[index]] = index
			}
			if(C$linkIndexB[index] > 0) {
				B$linkIndexC[C$linkIndexB[index]] = index
			}
			if(C$linkIndexD[index] > 0) {
				D$linkIndexC[C$linkIndexD[index]] = index
			}
		}

		# Sort by status
		D = D[order(D$status),]

		# Update possible links
		for(index in 1:nrow(D)) {
			if(D$linkIndexA[index] > 0) {
				A$linkIndexD[D$linkIndexA[index]] = index
			}
			if(D$linkIndexB[index] > 0) {
				B$linkIndexD[D$linkIndexB[index]] = index
			}
			if(D$linkIndexC[index] > 0) {
				C$linkIndexD[D$linkIndexC[index]] = index
			}
		}
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(B)) {
		cur = cur + B$relCount[index]
		B$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(C)) {
		cur = cur + C$relCount[index]
		C$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(D)) {
		cur = cur + D$relCount[index]
		D$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, status = "1A", accumRelCount = 0), A)
	dataB = rbind(data.frame(cloneCount = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, status = "1B", accumRelCount = 0), B)
	dataC = rbind(data.frame(cloneCount = 0, id = "START", factors = "C", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, status = "1C", accumRelCount = 0), C)
	dataD = rbind(data.frame(cloneCount = 0, id = "START", factors = "D", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, status = "1D", accumRelCount = 0), D)

	# Combine data sets
	data = rbind(dataA, dataB, dataC, dataD)

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}

	# Define link color ramp for each type of link
	linkRampAB = function(x) {
		ramp = colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampAC = function(x) {
		ramp = colorRamp(c(linkColors[3], linkColors[4]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampAD = function(x) {
		ramp = colorRamp(c(linkColors[5], linkColors[6]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampBC = function(x) {
		ramp = colorRamp(c(linkColors[7], linkColors[8]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampBD = function(x) {
		ramp = colorRamp(c(linkColors[9], linkColors[10]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}
	linkRampCD = function(x) {
		ramp = colorRamp(c(linkColors[11], linkColors[12]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
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

	# Add links
	for(index in 1:nrow(A)) {
		if(A$linkIndexB[index] > 0 && showLinks[1] == TRUE) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRampAB(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
		}
		if(A$linkIndexC[index] > 0 && showLinks[2] == TRUE) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "C", c(C$accumRelCount[A$linkIndexC[index]] - C$relCount[A$linkIndexC[index]], C$accumRelCount[A$linkIndexC[index]]), col = linkRampAC(A$relCount[index] + C$relCount[A$linkIndexC[index]]), border = NA)
		}
		if(A$linkIndexD[index] > 0 && showLinks[3] == TRUE) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "D", c(D$accumRelCount[A$linkIndexD[index]] - D$relCount[A$linkIndexD[index]], D$accumRelCount[A$linkIndexD[index]]), col = linkRampAD(A$relCount[index] + D$relCount[A$linkIndexD[index]]), border = NA)
		}
	}
	for(index in 1:nrow(B)) {
		if(B$linkIndexC[index] > 0 && showLinks[4] == TRUE) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "C", c(C$accumRelCount[B$linkIndexC[index]] - C$relCount[B$linkIndexC[index]], C$accumRelCount[B$linkIndexC[index]]), col = linkRampBC(B$relCount[index] + C$relCount[B$linkIndexC[index]]), border = NA)
		}
		if(B$linkIndexD[index] > 0 && showLinks[5] == TRUE) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "D", c(D$accumRelCount[B$linkIndexD[index]] - D$relCount[B$linkIndexD[index]], D$accumRelCount[B$linkIndexD[index]]), col = linkRampBD(B$relCount[index] + D$relCount[B$linkIndexD[index]]), border = NA)
		}
	}
	for(index in 1:nrow(C)) {
		if(C$linkIndexD[index] > 0 && showLinks[6] == TRUE) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "D", c(D$accumRelCount[C$linkIndexD[index]] - D$relCount[C$linkIndexD[index]], D$accumRelCount[C$linkIndexD[index]]), col = linkRampCD(C$relCount[index] + D$relCount[C$linkIndexD[index]]), border = NA)
		}
	}

	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos4.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			fileAlias4 = c(fileAliases[4]),
			cufoff = c(cutoff),
			sort = c(sort),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount),
			linkColorsAB1 = c(linkColors[1]),
			linkColorsAB2 = c(linkColors[2]),
			linkColorsAC1 = c(linkColors[3]),
			linkColorsAC2 = c(linkColors[4]),
			linkColorsAD1 = c(linkColors[5]),
			linkColorsAD2 = c(linkColors[6]),
			linkColorsBC1 = c(linkColors[7]),
			linkColorsBC2 = c(linkColors[8]),
			linkColorsBD1 = c(linkColors[9]),
			linkColorsBD2 = c(linkColors[10]),
			linkColorsCD1 = c(linkColors[11]),
			linkColorsCD2 = c(linkColors[12]),
			maxLinkSize = c(maxLinkSize),
			showLinksAB = c(showLinks[1]),
			showLinksAC = c(showLinks[2]),
			showLinksBC = c(showLinks[3]),
			showLinksAB = c(showLinks[4]),
			showLinksAC = c(showLinks[5]),
			showLinksBC = c(showLinks[6]))
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

circos6 = function(file1, file2, file3, file4, file5, file6, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080")) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 6) {
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
	if(length(linkColors) < 1) {
		print("Too few linkColors specified. Defaulting to red.")
		linkColors = rep(c("#FF000080"), each = 2)
	}
	if(length(linkColors) == 1) {
		linkColors = rep(linkColors, each = 2)
	}

	getStatus = function(data) {
		temp = data.frame(one = unlist(data[5:10]), two = c("A", "B", "C", "D", "E", "F"))
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
	file4Parts = unlist(strsplit(file4, "/"))
	file4Name = unlist(strsplit(file4Parts[length(file4Parts)], "[.]"))[1]
	file5Parts = unlist(strsplit(file5, "/"))
	file5Name = unlist(strsplit(file5Parts[length(file5Parts)], "[.]"))[1]
	file6Parts = unlist(strsplit(file6, "/"))
	file6Name = unlist(strsplit(file6Parts[length(file6Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name, file2Name, file3Name, file4Name, file5Name, file6Name)
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
	diagrammName = paste(file1Name, file2Name, file3Name, file4Name, file5Name, file6Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum
	A$linkIndexA = 1
	A$linkIndexB = 0
	A$linkIndexC = 0
	A$linkIndexD = 0
	A$linkIndexE = 0
	A$linkIndexF = 0

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

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = "B"
	sum = sum(B$cloneCount)
	B$relCount = B$cloneCount / sum
	B$linkIndexA = 0
	B$linkIndexB = 1
	B$linkIndexC = 0
	B$linkIndexD = 0
	B$linkIndexE = 0
	B$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		for(index in 1:nrow(B)) {
			cur = cur + B$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Read data from file3 and add additional columns
	C = readData(file3)
	C$factors = "C"
	sum = sum(C$cloneCount)
	C$relCount = C$cloneCount / sum
	C$linkIndexA = 0
	C$linkIndexB = 0
	C$linkIndexC = 1
	C$linkIndexD = 0
	C$linkIndexE = 0
	C$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(C)
		for(index in 1:nrow(C)) {
			cur = cur + C$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		C = C[1:cutoffIndex ,]
	}

	# Read data from file4 and add additional columns
	D = readData(file4)
	D$factors = "D"
	sum = sum(D$cloneCount)
	D$relCount = D$cloneCount / sum
	D$linkIndexA = 0
	D$linkIndexB = 0
	D$linkIndexC = 0
	D$linkIndexD = 1
	D$linkIndexE = 0
	D$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(D)
		for(index in 1:nrow(D)) {
			cur = cur + D$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		D = D[1:cutoffIndex ,]
	}

	# Read data from file5 and add additional columns
	E = readData(file5)
	E$factors = "E"
	sum = sum(E$cloneCount)
	E$relCount = E$cloneCount / sum
	E$linkIndexA = 0
	E$linkIndexB = 0
	E$linkIndexC = 0
	E$linkIndexD = 0
	E$linkIndexE = 1
	E$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(E)
		for(index in 1:nrow(E)) {
			cur = cur + E$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		E = E[1:cutoffIndex ,]
	}

	# Read data from file6 and add additional columns
	F = readData(file6)
	F$factors = "F"
	sum = sum(F$cloneCount)
	F$relCount = F$cloneCount / sum
	F$linkIndexA = 0
	F$linkIndexB = 0
	F$linkIndexC = 0
	F$linkIndexD = 0
	F$linkIndexE = 0
	F$linkIndexF = 1

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(F)
		for(index in 1:nrow(F)) {
			cur = cur + F$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		F = F[1:cutoffIndex ,]
	}

	# Compute links between A and B
	maxLinkSize = 0
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + B$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and C
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(C)) {
			if(A$id[index_1] == C$id[index_2]){
				A$linkIndexC[index_1] = index_2
				C$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and D
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(D)) {
			if(A$id[index_1] == D$id[index_2]){
				A$linkIndexD[index_1] = index_2
				D$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and E
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(E)) {
			if(A$id[index_1] == E$id[index_2]){
				A$linkIndexE[index_1] = index_2
				E$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and F
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(F)) {
			if(A$id[index_1] == F$id[index_2]){
				A$linkIndexF[index_1] = index_2
				F$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and C
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(C)) {
			if(B$id[index_1] == C$id[index_2]){
				B$linkIndexC[index_1] = index_2
				C$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and D
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(D)) {
			if(B$id[index_1] == D$id[index_2]){
				B$linkIndexD[index_1] = index_2
				D$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and E
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(E)) {
			if(B$id[index_1] == E$id[index_2]){
				B$linkIndexE[index_1] = index_2
				E$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and F
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(F)) {
			if(B$id[index_1] == F$id[index_2]){
				B$linkIndexF[index_1] = index_2
				F$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and D
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(D)) {
			if(C$id[index_1] == D$id[index_2]){
				C$linkIndexD[index_1] = index_2
				D$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and E
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(E)) {
			if(C$id[index_1] == E$id[index_2]){
				C$linkIndexE[index_1] = index_2
				E$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and F
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(F)) {
			if(C$id[index_1] == F$id[index_2]){
				C$linkIndexF[index_1] = index_2
				F$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between D and E
	for(index_1 in 1:nrow(D)) {
		for(index_2 in 1:nrow(E)) {
			if(D$id[index_1] == E$id[index_2]){
				D$linkIndexE[index_1] = index_2
				E$linkIndexD[index_2] = index_1
				linkSize = D$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between D and F
	for(index_1 in 1:nrow(D)) {
		for(index_2 in 1:nrow(F)) {
			if(D$id[index_1] == F$id[index_2]){
				D$linkIndexF[index_1] = index_2
				F$linkIndexD[index_2] = index_1
				linkSize = D$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between E and F
	for(index_1 in 1:nrow(E)) {
		for(index_2 in 1:nrow(F)) {
			if(E$id[index_1] == F$id[index_2]){
				E$linkIndexF[index_1] = index_2
				F$linkIndexE[index_2] = index_1
				linkSize = E$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
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

	# Set status
	for(index in 1:nrow(D)) {
		D$status[index] = getStatus(D[index,])
	}

	# Set status
	for(index in 1:nrow(E)) {
		E$status[index] = getStatus(E[index,])
	}

	# Set status
	for(index in 1:nrow(F)) {
		F$status[index] = getStatus(F[index,])
	}

	if(sort == TRUE) {
		# Sort by status
		A = A[order(A$status),]

		# Update possible links
		for(index in 1:nrow(A)) {
			if(A$linkIndexB[index] > 0) {
				B$linkIndexA[A$linkIndexB[index]] = index
			}
			if(A$linkIndexC[index] > 0) {
				C$linkIndexA[A$linkIndexC[index]] = index
			}
			if(A$linkIndexD[index] > 0) {
				D$linkIndexA[A$linkIndexD[index]] = index
			}
			if(A$linkIndexE[index] > 0) {
				E$linkIndexA[A$linkIndexE[index]] = index
			}
			if(A$linkIndexF[index] > 0) {
				F$linkIndexA[A$linkIndexF[index]] = index
			}
		}

		# Sort by status
		B = B[order(B$status),]

		# Update possible links
		for(index in 1:nrow(B)) {
			if(B$linkIndexA[index] > 0) {
				A$linkIndexB[B$linkIndexA[index]] = index
			}
			if(B$linkIndexC[index] > 0) {
				C$linkIndexB[B$linkIndexC[index]] = index
			}
			if(B$linkIndexD[index] > 0) {
				D$linkIndexB[B$linkIndexD[index]] = index
			}
			if(B$linkIndexE[index] > 0) {
				E$linkIndexB[B$linkIndexE[index]] = index
			}
			if(B$linkIndexF[index] > 0) {
				F$linkIndexB[B$linkIndexF[index]] = index
			}
		}

		# Sort by status
		C = C[order(C$status),]

		# Update possible links
		for(index in 1:nrow(C)) {
			if(C$linkIndexA[index] > 0) {
				A$linkIndexC[C$linkIndexA[index]] = index
			}
			if(C$linkIndexB[index] > 0) {
				B$linkIndexC[C$linkIndexB[index]] = index
			}
			if(C$linkIndexD[index] > 0) {
				D$linkIndexC[C$linkIndexD[index]] = index
			}
			if(C$linkIndexE[index] > 0) {
				E$linkIndexC[C$linkIndexE[index]] = index
			}
			if(C$linkIndexF[index] > 0) {
				F$linkIndexC[C$linkIndexF[index]] = index
			}
		}

		# Sort by status
		D = D[order(D$status),]

		# Update possible links
		for(index in 1:nrow(D)) {
			if(D$linkIndexA[index] > 0) {
				A$linkIndexD[D$linkIndexA[index]] = index
			}
			if(D$linkIndexB[index] > 0) {
				B$linkIndexD[D$linkIndexB[index]] = index
			}
			if(D$linkIndexC[index] > 0) {
				C$linkIndexD[D$linkIndexC[index]] = index
			}
			if(D$linkIndexE[index] > 0) {
				E$linkIndexD[D$linkIndexE[index]] = index
			}
			if(D$linkIndexF[index] > 0) {
				F$linkIndexD[D$linkIndexF[index]] = index
			}
		}

		# Sort by status
		E = E[order(E$status),]

		# Update possible links
		for(index in 1:nrow(E)) {
			if(E$linkIndexA[index] > 0) {
				A$linkIndexE[E$linkIndexA[index]] = index
			}
			if(E$linkIndexB[index] > 0) {
				B$linkIndexE[E$linkIndexB[index]] = index
			}
			if(E$linkIndexC[index] > 0) {
				C$linkIndexE[E$linkIndexC[index]] = index
			}
			if(E$linkIndexD[index] > 0) {
				D$linkIndexE[E$linkIndexD[index]] = index
			}
			if(E$linkIndexF[index] > 0) {
				F$linkIndexE[E$linkIndexF[index]] = index
			}
		}

		# Sort by status
		F = F[order(F$status),]

		# Update possible links
		for(index in 1:nrow(F)) {
			if(F$linkIndexA[index] > 0) {
				A$linkIndexF[F$linkIndexA[index]] = index
			}
			if(F$linkIndexB[index] > 0) {
				B$linkIndexF[F$linkIndexB[index]] = index
			}
			if(F$linkIndexC[index] > 0) {
				C$linkIndexF[F$linkIndexC[index]] = index
			}
			if(F$linkIndexD[index] > 0) {
				D$linkIndexF[F$linkIndexD[index]] = index
			}
			if(F$linkIndexE[index] > 0) {
				E$linkIndexF[F$linkIndexE[index]] = index
			}
		}
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(B)) {
		cur = cur + B$relCount[index]
		B$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(C)) {
		cur = cur + C$relCount[index]
		C$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(D)) {
		cur = cur + D$relCount[index]
		D$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(E)) {
		cur = cur + E$relCount[index]
		E$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(F)) {
		cur = cur + F$relCount[index]
		F$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1A", accumRelCount = 0), A)
	dataB = rbind(data.frame(cloneCount = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1B", accumRelCount = 0), B)
	dataC = rbind(data.frame(cloneCount = 0, id = "START", factors = "C", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1C", accumRelCount = 0), C)
	dataD = rbind(data.frame(cloneCount = 0, id = "START", factors = "D", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1D", accumRelCount = 0), D)
	dataE = rbind(data.frame(cloneCount = 0, id = "START", factors = "E", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1E", accumRelCount = 0), E)
	dataF = rbind(data.frame(cloneCount = 0, id = "START", factors = "F", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1F", accumRelCount = 0), F)

	# Combine data sets
	data = rbind(dataA, dataB, dataC, dataD, dataE, dataF)

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}

	# Define link color ramp
	linkRamp = function(x) {
		ramp = colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
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

	# Add links
	for(index in 1:nrow(A)) {
		if(A$linkIndexB[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRamp(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
		}
		if(A$linkIndexC[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "C", c(C$accumRelCount[A$linkIndexC[index]] - C$relCount[A$linkIndexC[index]], C$accumRelCount[A$linkIndexC[index]]), col = linkRamp(A$relCount[index] + C$relCount[A$linkIndexC[index]]), border = NA)
		}
		if(A$linkIndexD[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "D", c(D$accumRelCount[A$linkIndexD[index]] - D$relCount[A$linkIndexD[index]], D$accumRelCount[A$linkIndexD[index]]), col = linkRamp(A$relCount[index] + D$relCount[A$linkIndexD[index]]), border = NA)
		}
		if(A$linkIndexE[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "E", c(E$accumRelCount[A$linkIndexE[index]] - E$relCount[A$linkIndexE[index]], E$accumRelCount[A$linkIndexE[index]]), col = linkRamp(A$relCount[index] + E$relCount[A$linkIndexE[index]]), border = NA)
		}
		if(A$linkIndexF[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "F", c(F$accumRelCount[A$linkIndexF[index]] - F$relCount[A$linkIndexF[index]], F$accumRelCount[A$linkIndexF[index]]), col = linkRamp(A$relCount[index] + F$relCount[A$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(B)) {
		if(B$linkIndexC[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "C", c(C$accumRelCount[B$linkIndexC[index]] - C$relCount[B$linkIndexC[index]], C$accumRelCount[B$linkIndexC[index]]), col = linkRamp(B$relCount[index] + C$relCount[B$linkIndexC[index]]), border = NA)
		}
		if(B$linkIndexD[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "D", c(D$accumRelCount[B$linkIndexD[index]] - D$relCount[B$linkIndexD[index]], D$accumRelCount[B$linkIndexD[index]]), col = linkRamp(B$relCount[index] + D$relCount[B$linkIndexD[index]]), border = NA)
		}
		if(B$linkIndexE[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "E", c(E$accumRelCount[B$linkIndexE[index]] - E$relCount[B$linkIndexE[index]], E$accumRelCount[B$linkIndexE[index]]), col = linkRamp(B$relCount[index] + E$relCount[B$linkIndexE[index]]), border = NA)
		}
		if(B$linkIndexF[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "F", c(F$accumRelCount[B$linkIndexF[index]] - F$relCount[B$linkIndexF[index]], F$accumRelCount[B$linkIndexF[index]]), col = linkRamp(B$relCount[index] + F$relCount[B$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(C)) {
		if(C$linkIndexD[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "D", c(D$accumRelCount[C$linkIndexD[index]] - D$relCount[C$linkIndexD[index]], D$accumRelCount[C$linkIndexD[index]]), col = linkRamp(C$relCount[index] + D$relCount[C$linkIndexD[index]]), border = NA)
		}
		if(C$linkIndexE[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "E", c(E$accumRelCount[C$linkIndexE[index]] - E$relCount[C$linkIndexE[index]], E$accumRelCount[C$linkIndexE[index]]), col = linkRamp(C$relCount[index] + E$relCount[C$linkIndexE[index]]), border = NA)
		}
		if(C$linkIndexF[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "F", c(F$accumRelCount[C$linkIndexF[index]] - F$relCount[C$linkIndexF[index]], F$accumRelCount[C$linkIndexF[index]]), col = linkRamp(C$relCount[index] + F$relCount[C$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(D)) {
		if(D$linkIndexE[index] > 0) {
			circos.link("D", c(D$accumRelCount[index] - D$relCount[index], D$accumRelCount[index]), "E", c(E$accumRelCount[D$linkIndexE[index]] - E$relCount[D$linkIndexE[index]], E$accumRelCount[D$linkIndexE[index]]), col = linkRamp(D$relCount[index] + E$relCount[D$linkIndexE[index]]), border = NA)
		}
		if(D$linkIndexF[index] > 0) {
			circos.link("D", c(D$accumRelCount[index] - D$relCount[index], D$accumRelCount[index]), "F", c(F$accumRelCount[D$linkIndexF[index]] - F$relCount[D$linkIndexF[index]], F$accumRelCount[D$linkIndexF[index]]), col = linkRamp(D$relCount[index] + F$relCount[D$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(E)) {
		if(E$linkIndexF[index] > 0) {
			circos.link("E", c(E$accumRelCount[index] - E$relCount[index], E$accumRelCount[index]), "F", c(F$accumRelCount[E$linkIndexF[index]] - F$relCount[E$linkIndexF[index]], F$accumRelCount[E$linkIndexF[index]]), col = linkRamp(E$relCount[index] + F$relCount[E$linkIndexF[index]]), border = NA)
		}
	}
	
	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos6.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			fileAlias4 = c(fileAliases[4]),
			fileAlias5 = c(fileAliases[5]),
			fileAlias6 = c(fileAliases[6]),
			cufoff = c(cutoff),
			sort = c(sort),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount),
			linkColors1 = c(linkColors[1]),
			linkColors2 = c(linkColors[2]),
			maxLinkSize = c(maxLinkSize))
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


circos6 = function(file1, file2, file3, file4, file5, file6, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080")) {
	# Catch invalid argument values
	if(!is.null(fileAliases) && length(fileAliases) < 6) {
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
	if(length(linkColors) < 1) {
		print("Too few linkColors specified. Defaulting to red.")
		linkColors = rep(c("#FF000080"), each = 2)
	}
	if(length(linkColors) == 1) {
		linkColors = rep(linkColors, each = 2)
	}

	getStatus = function(data) {
		temp = data.frame(one = unlist(data[5:10]), two = c("A", "B", "C", "D", "E", "F"))
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
	file4Parts = unlist(strsplit(file4, "/"))
	file4Name = unlist(strsplit(file4Parts[length(file4Parts)], "[.]"))[1]
	file5Parts = unlist(strsplit(file5, "/"))
	file5Name = unlist(strsplit(file5Parts[length(file5Parts)], "[.]"))[1]
	file6Parts = unlist(strsplit(file6, "/"))
	file6Name = unlist(strsplit(file6Parts[length(file6Parts)], "[.]"))[1]

	# Set file aliases if not provided
	if(is.null(fileAliases)) {
		fileAliases = c(file1Name, file2Name, file3Name, file4Name, file5Name, file6Name)
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
	diagrammName = paste(file1Name, file2Name, file3Name, file4Name, file5Name, file6Name, sep = "_")
	diagrammFileName = paste(saveFolder, diagrammName, "_circos.svg", sep = "")

	# Read data from file1 and add additional columns
	A = readData(file1)
	A$factors = "A"
	sum = sum(A$cloneCount)
	A$relCount = A$cloneCount / sum
	A$linkIndexA = 1
	A$linkIndexB = 0
	A$linkIndexC = 0
	A$linkIndexD = 0
	A$linkIndexE = 0
	A$linkIndexF = 0

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

	# Read data from file2 and add additional columns
	B = readData(file2)
	B$factors = "B"
	sum = sum(B$cloneCount)
	B$relCount = B$cloneCount / sum
	B$linkIndexA = 0
	B$linkIndexB = 1
	B$linkIndexC = 0
	B$linkIndexD = 0
	B$linkIndexE = 0
	B$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(B)
		for(index in 1:nrow(B)) {
			cur = cur + B$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		B = B[1:cutoffIndex ,]
	}

	# Read data from file3 and add additional columns
	C = readData(file3)
	C$factors = "C"
	sum = sum(C$cloneCount)
	C$relCount = C$cloneCount / sum
	C$linkIndexA = 0
	C$linkIndexB = 0
	C$linkIndexC = 1
	C$linkIndexD = 0
	C$linkIndexE = 0
	C$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(C)
		for(index in 1:nrow(C)) {
			cur = cur + C$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		C = C[1:cutoffIndex ,]
	}

	# Read data from file4 and add additional columns
	D = readData(file4)
	D$factors = "D"
	sum = sum(D$cloneCount)
	D$relCount = D$cloneCount / sum
	D$linkIndexA = 0
	D$linkIndexB = 0
	D$linkIndexC = 0
	D$linkIndexD = 1
	D$linkIndexE = 0
	D$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(D)
		for(index in 1:nrow(D)) {
			cur = cur + D$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		D = D[1:cutoffIndex ,]
	}

	# Read data from file5 and add additional columns
	E = readData(file5)
	E$factors = "E"
	sum = sum(E$cloneCount)
	E$relCount = E$cloneCount / sum
	E$linkIndexA = 0
	E$linkIndexB = 0
	E$linkIndexC = 0
	E$linkIndexD = 0
	E$linkIndexE = 1
	E$linkIndexF = 0

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(E)
		for(index in 1:nrow(E)) {
			cur = cur + E$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		E = E[1:cutoffIndex ,]
	}

	# Read data from file6 and add additional columns
	F = readData(file6)
	F$factors = "F"
	sum = sum(F$cloneCount)
	F$relCount = F$cloneCount / sum
	F$linkIndexA = 0
	F$linkIndexB = 0
	F$linkIndexC = 0
	F$linkIndexD = 0
	F$linkIndexE = 0
	F$linkIndexF = 1

	# Perform cutoff
	if(cutoff < 1) {
		cur = 0
		cutoffIndex = nrow(F)
		for(index in 1:nrow(F)) {
			cur = cur + F$relCount[index]
			if(cur > cutoff) {
				cutoffIndex = index
				break
			}
		}
		F = F[1:cutoffIndex ,]
	}

	# Compute links between A and B
	maxLinkSize = 0
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(B)) {
			if(A$id[index_1] == B$id[index_2]){
				A$linkIndexB[index_1] = index_2
				B$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + B$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and C
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(C)) {
			if(A$id[index_1] == C$id[index_2]){
				A$linkIndexC[index_1] = index_2
				C$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and D
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(D)) {
			if(A$id[index_1] == D$id[index_2]){
				A$linkIndexD[index_1] = index_2
				D$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and E
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(E)) {
			if(A$id[index_1] == E$id[index_2]){
				A$linkIndexE[index_1] = index_2
				E$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between A and F
	for(index_1 in 1:nrow(A)) {
		for(index_2 in 1:nrow(F)) {
			if(A$id[index_1] == F$id[index_2]){
				A$linkIndexF[index_1] = index_2
				F$linkIndexA[index_2] = index_1
				linkSize = A$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and C
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(C)) {
			if(B$id[index_1] == C$id[index_2]){
				B$linkIndexC[index_1] = index_2
				C$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + C$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and D
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(D)) {
			if(B$id[index_1] == D$id[index_2]){
				B$linkIndexD[index_1] = index_2
				D$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and E
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(E)) {
			if(B$id[index_1] == E$id[index_2]){
				B$linkIndexE[index_1] = index_2
				E$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between B and F
	for(index_1 in 1:nrow(B)) {
		for(index_2 in 1:nrow(F)) {
			if(B$id[index_1] == F$id[index_2]){
				B$linkIndexF[index_1] = index_2
				F$linkIndexB[index_2] = index_1
				linkSize = B$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and D
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(D)) {
			if(C$id[index_1] == D$id[index_2]){
				C$linkIndexD[index_1] = index_2
				D$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + D$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and E
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(E)) {
			if(C$id[index_1] == E$id[index_2]){
				C$linkIndexE[index_1] = index_2
				E$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between C and F
	for(index_1 in 1:nrow(C)) {
		for(index_2 in 1:nrow(F)) {
			if(C$id[index_1] == F$id[index_2]){
				C$linkIndexF[index_1] = index_2
				F$linkIndexC[index_2] = index_1
				linkSize = C$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between D and E
	for(index_1 in 1:nrow(D)) {
		for(index_2 in 1:nrow(E)) {
			if(D$id[index_1] == E$id[index_2]){
				D$linkIndexE[index_1] = index_2
				E$linkIndexD[index_2] = index_1
				linkSize = D$relCount[index_1] + E$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between D and F
	for(index_1 in 1:nrow(D)) {
		for(index_2 in 1:nrow(F)) {
			if(D$id[index_1] == F$id[index_2]){
				D$linkIndexF[index_1] = index_2
				F$linkIndexD[index_2] = index_1
				linkSize = D$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
			}
		}
	}

	# Compute links between E and F
	for(index_1 in 1:nrow(E)) {
		for(index_2 in 1:nrow(F)) {
			if(E$id[index_1] == F$id[index_2]){
				E$linkIndexF[index_1] = index_2
				F$linkIndexE[index_2] = index_1
				linkSize = E$relCount[index_1] + F$relCount[index_2]
				maxLinkSize = ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
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

	# Set status
	for(index in 1:nrow(D)) {
		D$status[index] = getStatus(D[index,])
	}

	# Set status
	for(index in 1:nrow(E)) {
		E$status[index] = getStatus(E[index,])
	}

	# Set status
	for(index in 1:nrow(F)) {
		F$status[index] = getStatus(F[index,])
	}

	if(sort == TRUE) {
		# Sort by status
		A = A[order(A$status),]

		# Update possible links
		for(index in 1:nrow(A)) {
			if(A$linkIndexB[index] > 0) {
				B$linkIndexA[A$linkIndexB[index]] = index
			}
			if(A$linkIndexC[index] > 0) {
				C$linkIndexA[A$linkIndexC[index]] = index
			}
			if(A$linkIndexD[index] > 0) {
				D$linkIndexA[A$linkIndexD[index]] = index
			}
			if(A$linkIndexE[index] > 0) {
				E$linkIndexA[A$linkIndexE[index]] = index
			}
			if(A$linkIndexF[index] > 0) {
				F$linkIndexA[A$linkIndexF[index]] = index
			}
		}

		# Sort by status
		B = B[order(B$status),]

		# Update possible links
		for(index in 1:nrow(B)) {
			if(B$linkIndexA[index] > 0) {
				A$linkIndexB[B$linkIndexA[index]] = index
			}
			if(B$linkIndexC[index] > 0) {
				C$linkIndexB[B$linkIndexC[index]] = index
			}
			if(B$linkIndexD[index] > 0) {
				D$linkIndexB[B$linkIndexD[index]] = index
			}
			if(B$linkIndexE[index] > 0) {
				E$linkIndexB[B$linkIndexE[index]] = index
			}
			if(B$linkIndexF[index] > 0) {
				F$linkIndexB[B$linkIndexF[index]] = index
			}
		}

		# Sort by status
		C = C[order(C$status),]

		# Update possible links
		for(index in 1:nrow(C)) {
			if(C$linkIndexA[index] > 0) {
				A$linkIndexC[C$linkIndexA[index]] = index
			}
			if(C$linkIndexB[index] > 0) {
				B$linkIndexC[C$linkIndexB[index]] = index
			}
			if(C$linkIndexD[index] > 0) {
				D$linkIndexC[C$linkIndexD[index]] = index
			}
			if(C$linkIndexE[index] > 0) {
				E$linkIndexC[C$linkIndexE[index]] = index
			}
			if(C$linkIndexF[index] > 0) {
				F$linkIndexC[C$linkIndexF[index]] = index
			}
		}

		# Sort by status
		D = D[order(D$status),]

		# Update possible links
		for(index in 1:nrow(D)) {
			if(D$linkIndexA[index] > 0) {
				A$linkIndexD[D$linkIndexA[index]] = index
			}
			if(D$linkIndexB[index] > 0) {
				B$linkIndexD[D$linkIndexB[index]] = index
			}
			if(D$linkIndexC[index] > 0) {
				C$linkIndexD[D$linkIndexC[index]] = index
			}
			if(D$linkIndexE[index] > 0) {
				E$linkIndexD[D$linkIndexE[index]] = index
			}
			if(D$linkIndexF[index] > 0) {
				F$linkIndexD[D$linkIndexF[index]] = index
			}
		}

		# Sort by status
		E = E[order(E$status),]

		# Update possible links
		for(index in 1:nrow(E)) {
			if(E$linkIndexA[index] > 0) {
				A$linkIndexE[E$linkIndexA[index]] = index
			}
			if(E$linkIndexB[index] > 0) {
				B$linkIndexE[E$linkIndexB[index]] = index
			}
			if(E$linkIndexC[index] > 0) {
				C$linkIndexE[E$linkIndexC[index]] = index
			}
			if(E$linkIndexD[index] > 0) {
				D$linkIndexE[E$linkIndexD[index]] = index
			}
			if(E$linkIndexF[index] > 0) {
				F$linkIndexE[E$linkIndexF[index]] = index
			}
		}

		# Sort by status
		F = F[order(F$status),]

		# Update possible links
		for(index in 1:nrow(F)) {
			if(F$linkIndexA[index] > 0) {
				A$linkIndexF[F$linkIndexA[index]] = index
			}
			if(F$linkIndexB[index] > 0) {
				B$linkIndexF[F$linkIndexB[index]] = index
			}
			if(F$linkIndexC[index] > 0) {
				C$linkIndexF[F$linkIndexC[index]] = index
			}
			if(F$linkIndexD[index] > 0) {
				D$linkIndexF[F$linkIndexD[index]] = index
			}
			if(F$linkIndexE[index] > 0) {
				E$linkIndexF[F$linkIndexE[index]] = index
			}
		}
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(A)) {
		cur = cur + A$relCount[index]
		A$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(B)) {
		cur = cur + B$relCount[index]
		B$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(C)) {
		cur = cur + C$relCount[index]
		C$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(D)) {
		cur = cur + D$relCount[index]
		D$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(E)) {
		cur = cur + E$relCount[index]
		E$accumRelCount[index] = cur
	}

	# Accumulate relative counts
	cur = 0
	for(index in 1:nrow(F)) {
		cur = cur + F$relCount[index]
		F$accumRelCount[index] = cur
	}

	# Add a starting row to each data set
	dataA = rbind(data.frame(cloneCount = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1A", accumRelCount = 0), A)
	dataB = rbind(data.frame(cloneCount = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1B", accumRelCount = 0), B)
	dataC = rbind(data.frame(cloneCount = 0, id = "START", factors = "C", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1C", accumRelCount = 0), C)
	dataD = rbind(data.frame(cloneCount = 0, id = "START", factors = "D", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1D", accumRelCount = 0), D)
	dataE = rbind(data.frame(cloneCount = 0, id = "START", factors = "E", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1E", accumRelCount = 0), E)
	dataF = rbind(data.frame(cloneCount = 0, id = "START", factors = "F", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, linkIndexD = 0, linkIndexE = 0, linkIndexF = 0, status = "1F", accumRelCount = 0), F)

	# Combine data sets
	data = rbind(dataA, dataB, dataC, dataD, dataE, dataF)

	# Compute size of largest clone
	maxCount = max(data$relCount)

	# Define count color ramp
	countRamp = function(x) {
		ramp = colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
		color = ramp(x/maxCount)
		rgb(color, alpha = color[4], maxColorValue =  255)
	}

	# Define link color ramp
	linkRamp = function(x) {
		ramp = colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
		color = ramp(x/maxLinkSize)
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

	# Add links
	for(index in 1:nrow(A)) {
		if(A$linkIndexB[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRamp(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
		}
		if(A$linkIndexC[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "C", c(C$accumRelCount[A$linkIndexC[index]] - C$relCount[A$linkIndexC[index]], C$accumRelCount[A$linkIndexC[index]]), col = linkRamp(A$relCount[index] + C$relCount[A$linkIndexC[index]]), border = NA)
		}
		if(A$linkIndexD[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "D", c(D$accumRelCount[A$linkIndexD[index]] - D$relCount[A$linkIndexD[index]], D$accumRelCount[A$linkIndexD[index]]), col = linkRamp(A$relCount[index] + D$relCount[A$linkIndexD[index]]), border = NA)
		}
		if(A$linkIndexE[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "E", c(E$accumRelCount[A$linkIndexE[index]] - E$relCount[A$linkIndexE[index]], E$accumRelCount[A$linkIndexE[index]]), col = linkRamp(A$relCount[index] + E$relCount[A$linkIndexE[index]]), border = NA)
		}
		if(A$linkIndexF[index] > 0) {
			circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "F", c(F$accumRelCount[A$linkIndexF[index]] - F$relCount[A$linkIndexF[index]], F$accumRelCount[A$linkIndexF[index]]), col = linkRamp(A$relCount[index] + F$relCount[A$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(B)) {
		if(B$linkIndexC[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "C", c(C$accumRelCount[B$linkIndexC[index]] - C$relCount[B$linkIndexC[index]], C$accumRelCount[B$linkIndexC[index]]), col = linkRamp(B$relCount[index] + C$relCount[B$linkIndexC[index]]), border = NA)
		}
		if(B$linkIndexD[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "D", c(D$accumRelCount[B$linkIndexD[index]] - D$relCount[B$linkIndexD[index]], D$accumRelCount[B$linkIndexD[index]]), col = linkRamp(B$relCount[index] + D$relCount[B$linkIndexD[index]]), border = NA)
		}
		if(B$linkIndexE[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "E", c(E$accumRelCount[B$linkIndexE[index]] - E$relCount[B$linkIndexE[index]], E$accumRelCount[B$linkIndexE[index]]), col = linkRamp(B$relCount[index] + E$relCount[B$linkIndexE[index]]), border = NA)
		}
		if(B$linkIndexF[index] > 0) {
			circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "F", c(F$accumRelCount[B$linkIndexF[index]] - F$relCount[B$linkIndexF[index]], F$accumRelCount[B$linkIndexF[index]]), col = linkRamp(B$relCount[index] + F$relCount[B$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(C)) {
		if(C$linkIndexD[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "D", c(D$accumRelCount[C$linkIndexD[index]] - D$relCount[C$linkIndexD[index]], D$accumRelCount[C$linkIndexD[index]]), col = linkRamp(C$relCount[index] + D$relCount[C$linkIndexD[index]]), border = NA)
		}
		if(C$linkIndexE[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "E", c(E$accumRelCount[C$linkIndexE[index]] - E$relCount[C$linkIndexE[index]], E$accumRelCount[C$linkIndexE[index]]), col = linkRamp(C$relCount[index] + E$relCount[C$linkIndexE[index]]), border = NA)
		}
		if(C$linkIndexF[index] > 0) {
			circos.link("C", c(C$accumRelCount[index] - C$relCount[index], C$accumRelCount[index]), "F", c(F$accumRelCount[C$linkIndexF[index]] - F$relCount[C$linkIndexF[index]], F$accumRelCount[C$linkIndexF[index]]), col = linkRamp(C$relCount[index] + F$relCount[C$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(D)) {
		if(D$linkIndexE[index] > 0) {
			circos.link("D", c(D$accumRelCount[index] - D$relCount[index], D$accumRelCount[index]), "E", c(E$accumRelCount[D$linkIndexE[index]] - E$relCount[D$linkIndexE[index]], E$accumRelCount[D$linkIndexE[index]]), col = linkRamp(D$relCount[index] + E$relCount[D$linkIndexE[index]]), border = NA)
		}
		if(D$linkIndexF[index] > 0) {
			circos.link("D", c(D$accumRelCount[index] - D$relCount[index], D$accumRelCount[index]), "F", c(F$accumRelCount[D$linkIndexF[index]] - F$relCount[D$linkIndexF[index]], F$accumRelCount[D$linkIndexF[index]]), col = linkRamp(D$relCount[index] + F$relCount[D$linkIndexF[index]]), border = NA)
		}
	}
	for(index in 1:nrow(E)) {
		if(E$linkIndexF[index] > 0) {
			circos.link("E", c(E$accumRelCount[index] - E$relCount[index], E$accumRelCount[index]), "F", c(F$accumRelCount[E$linkIndexF[index]] - F$relCount[E$linkIndexF[index]], F$accumRelCount[E$linkIndexF[index]]), col = linkRamp(E$relCount[index] + F$relCount[E$linkIndexF[index]]), border = NA)
		}
	}
	
	# Close diagramm and output file
	circos.clear()
	dev.off()

	# Store parameters for later reference
	parameterFileName = paste(saveFolder, "circos6.csv", sep = "")
	parameter = data.frame(diagrammName = c(diagrammName),
			fileAlias1 = c(fileAliases[1]),
			fileAlias2 = c(fileAliases[2]),
			fileAlias3 = c(fileAliases[3]),
			fileAlias4 = c(fileAliases[4]),
			fileAlias5 = c(fileAliases[5]),
			fileAlias6 = c(fileAliases[6]),
			cufoff = c(cutoff),
			sort = c(sort),
			countColors1 = c(countColors[1]),
			countColors2 = c(countColors[2]),
			maxCount = c(maxCount),
			linkColors1 = c(linkColors[1]),
			linkColors2 = c(linkColors[2]),
			maxLinkSize = c(maxLinkSize))
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
