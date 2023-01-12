library(circlize)
circos3 <- function(file1, file2, file3, fileAliases = NULL, saveFolder = NULL, cutoff = 1.0, sort = FALSE, countColors = c("#FFFFFFFF", "#0000FFFF"), linkColors = c("#FF000080", "#FF000080", "#FF000080"), showLinks = c(TRUE, TRUE, TRUE)) {
    # Catch invalid argument values
    if (!is.null(fileAliases) && length(fileAliases) < 3) {
        print("To few file aliases specified. Defaulting to file names.")
        fileAliases <- NULL
    }
    if (cutoff < 0.0) {
        pirnt("Cutoff value is below zero. Defaulting to 1.0.")
        cutoff <- 1.0
    }
    if (cutoff > 1.0) {
        print("Cutoff value is above 1.0. Defaulting to 1.0.")
        cutoff <- 1.0
    }
    if (length(countColors) < 2) {
        print("Too few countColors specified. Defaulting to blue.")
        countColors <- c("#FFFFFFFF", "#0000FFFF")
    }
    if (length(linkColors) < 3) {
        print("Too few linkColors specified. Defaulting to red.")
        linkColors <- rep(c("#FF000080", "#FF000080", "#FF000080"), each = 2)
    }
    if (length(linkColors) == 3) {
        linkColors <- rep(linkColors, each = 2)
    }
    if (length(linkColors) > 3 && length(linkColors) < 6) {
        print("Too few linkColors for link color ramp specified. Defaulting to red.")
        linkColors <- rep(c("#FF000080", "#FF000080", "#FF000080"), each = 2)
    }
    if (length(showLinks) < 3) {
        print("Too few values in showLinks specified. Defaulting to showing all links.")
        showLinks <- c(TRUE, TRUE, TRUE)
    }

    # Helper functions
    readData <- function(file) {
        data <- read.table(file, header = TRUE, sep = "\t", colClasses = c("integer", "character", "character", "character", "character", "numeric"))
        data <- cbind(data, data.frame(id = paste(data$CDR3.nucleotide.sequence, data$V.segments, data$J.segments, sep = "")))
        subset(data, select = -c(CDR3.nucleotide.sequence, CDR3.amino.acid.sequence, V.segments, J.segments, Percentage))
    }
    getStatus <- function(data) {
        temp <- data.frame(one = unlist(data[5:7]), two = c("A", "B", "C"))
        statusLetters <- paste(temp[temp$one > 0, "two"], collapse = "")
        paste(nchar(statusLetters), statusLetters, sep = "")
    }

    # Extract file names
    file1Parts <- unlist(strsplit(file1, "/"))
    file1Name <- unlist(strsplit(file1Parts[length(file1Parts)], "[.]"))[1]
    file2Parts <- unlist(strsplit(file2, "/"))
    file2Name <- unlist(strsplit(file2Parts[length(file2Parts)], "[.]"))[1]
    file3Parts <- unlist(strsplit(file3, "/"))
    file3Name <- unlist(strsplit(file3Parts[length(file3Parts)], "[.]"))[1]

    # Set file aliases if not provided
    if (is.null(fileAliases)) {
        fileAliases <- c(file1Name, file2Name, file3Name)
    }

    # Create save folder name if not provided
    if (is.null(saveFolder)) {
        if (length(file1Parts) > 1) {
            saveFolder <- paste(file1Parts[1:length(file1Parts) - 1], collapse = "/")
            saveFolder <- paste(saveFolder, "/", sep = "")
        } else {
            saveFolder <- ""
        }
    }

    # Create diagramm name
    diagrammName <- paste(file1Name, file2Name, file3Name, sep = "_")
    diagrammFileName <- paste(saveFolder, diagrammName, "_circos.svg", sep = "")

    # Read data from file1 and add additional columns
    A <- readData(file1)
    A$factors <- "A"
    sum <- sum(A$Count)
    A$relCount <- A$Count / sum
    A$linkIndexA <- 1
    A$linkIndexB <- 0
    A$linkIndexC <- 0

    # Perform cutoff
    if (cutoff < 1) {
        cur <- 0
        cutoffIndex <- nrow(A)
        for (index in 1:nrow(A)) {
            cur <- cur + A$relCount[index]
            if (cur > cutoff) {
                cutoffIndex <- index
                break
            }
        }
        A <- A[1:cutoffIndex, ]
    }

    # Read data from file2 and add additional columns
    B <- readData(file2)
    B$factors <- "B"
    sum <- sum(B$Count)
    B$relCount <- B$Count / sum
    B$linkIndexA <- 0
    B$linkIndexB <- 1
    B$linkIndexC <- 0

    # Perform cutoff
    if (cutoff < 1) {
        cur <- 0
        cutoffIndex <- nrow(B)
        for (index in 1:nrow(B)) {
            cur <- cur + B$relCount[index]
            if (cur > cutoff) {
                cutoffIndex <- index
                break
            }
        }
        B <- B[1:cutoffIndex, ]
    }

    # Read data from file3 and add additional columns
    C <- readData(file3)
    C$factors <- "C"
    sum <- sum(C$Count)
    C$relCount <- C$Count / sum
    C$linkIndexA <- 0
    C$linkIndexB <- 0
    C$linkIndexC <- 1

    # Perform cutoff
    if (cutoff < 1) {
        cur <- 0
        cutoffIndex <- nrow(C)
        for (index in 1:nrow(C)) {
            cur <- cur + C$relCount[index]
            if (cur > cutoff) {
                cutoffIndex <- index
                break
            }
        }
        C <- C[1:cutoffIndex, ]
    }

    # Compute links between A and B
    maxLinkSize <- 0
    for (index_1 in 1:nrow(A)) {
        for (index_2 in 1:nrow(B)) {
            if (A$id[index_1] == B$id[index_2]) {
                A$linkIndexB[index_1] <- index_2
                B$linkIndexA[index_2] <- index_1
                linkSize <- A$relCount[index_1] + B$relCount[index_2]
                maxLinkSize <- ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
            }
        }
    }

    # Compute links between A and C
    for (index_1 in 1:nrow(A)) {
        for (index_2 in 1:nrow(C)) {
            if (A$id[index_1] == C$id[index_2]) {
                A$linkIndexC[index_1] <- index_2
                C$linkIndexA[index_2] <- index_1
                linkSize <- A$relCount[index_1] + C$relCount[index_2]
                maxLinkSize <- ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
            }
        }
    }

    # Compute links between B and C
    for (index_1 in 1:nrow(B)) {
        for (index_2 in 1:nrow(C)) {
            if (B$id[index_1] == C$id[index_2]) {
                B$linkIndexC[index_1] <- index_2
                C$linkIndexB[index_2] <- index_1
                linkSize <- B$relCount[index_1] + C$relCount[index_2]
                maxLinkSize <- ifelse(maxLinkSize < linkSize, linkSize, maxLinkSize)
            }
        }
    }

    # Set status
    for (index in 1:nrow(A)) {
        A$status[index] <- getStatus(A[index, ])
    }

    # Set status
    for (index in 1:nrow(B)) {
        B$status[index] <- getStatus(B[index, ])
    }

    # Set status
    for (index in 1:nrow(C)) {
        C$status[index] <- getStatus(C[index, ])
    }

    if (sort == TRUE) {
        # Sort by status
        A <- A[order(A$status), ]

        # Update possible links
        for (index in 1:nrow(A)) {
            if (A$linkIndexB[index] > 0) {
                B$linkIndexA[A$linkIndexB[index]] <- index
            }
            if (A$linkIndexC[index] > 0) {
                C$linkIndexA[A$linkIndexC[index]] <- index
            }
        }

        # Sort by status
        B <- B[order(B$status), ]

        # Update possible links
        for (index in 1:nrow(B)) {
            if (B$linkIndexA[index] > 0) {
                A$linkIndexB[B$linkIndexA[index]] <- index
            }
            if (B$linkIndexC[index] > 0) {
                C$linkIndexB[B$linkIndexC[index]] <- index
            }
        }

        # Sort by status
        C <- C[order(C$status), ]

        # Update possible links
        for (index in 1:nrow(C)) {
            if (C$linkIndexA[index] > 0) {
                A$linkIndexC[C$linkIndexA[index]] <- index
            }
            if (C$linkIndexB[index] > 0) {
                B$linkIndexC[C$linkIndexB[index]] <- index
            }
        }
    }

    # Accumulate relative counts
    cur <- 0
    for (index in 1:nrow(A)) {
        cur <- cur + A$relCount[index]
        A$accumRelCount[index] <- cur
    }

    # Accumulate relative counts
    cur <- 0
    for (index in 1:nrow(B)) {
        cur <- cur + B$relCount[index]
        B$accumRelCount[index] <- cur
    }

    # Accumulate relative counts
    cur <- 0
    for (index in 1:nrow(C)) {
        cur <- cur + C$relCount[index]
        C$accumRelCount[index] <- cur
    }

    # Add a starting row to each data set
    dataA <- rbind(data.frame(Count = 0, id = "START", factors = "A", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1A", accumRelCount = 0), A)
    dataB <- rbind(data.frame(Count = 0, id = "START", factors = "B", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1B", accumRelCount = 0), B)
    dataC <- rbind(data.frame(Count = 0, id = "START", factors = "C", relCount = 0, linkIndexA = 0, linkIndexB = 0, linkIndexC = 0, status = "1C", accumRelCount = 0), C)

    # Combine data sets
    data <- rbind(dataA, dataB, dataC)

    # Compute size of largest clone
    maxCount <- max(data$relCount)

    # Define count color ramp
    countRamp <- function(x) {
        ramp <- colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
        color <- ramp(x / maxCount)
        rgb(color, alpha = color[4], maxColorValue = 255)
    }

    # Define link color ramp for each type of link
    linkRampAB <- function(x) {
        ramp <- colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
        color <- ramp(x / maxLinkSize)
        rgb(color, alpha = color[4], maxColorValue = 255)
    }
    linkRampAC <- function(x) {
        ramp <- colorRamp(c(linkColors[3], linkColors[4]), alpha = TRUE)
        color <- ramp(x / maxLinkSize)
        rgb(color, alpha = color[4], maxColorValue = 255)
    }
    linkRampBC <- function(x) {
        ramp <- colorRamp(c(linkColors[5], linkColors[6]), alpha = TRUE)
        color <- ramp(x / maxLinkSize)
        rgb(color, alpha = color[4], maxColorValue = 255)
    }

    # Open output file, set parameters and initialize diagramm
    svg(diagrammFileName)
    circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
    circos.initialize(factors = data$factors, x = data$accumRelCount)

    # Add ring with labels, ticks and colored boxes
    circos.track(
        factors = data$factors, x = data$relCount, y = data$accumRelCount, bg.border = NA,
        panel.fun = function(x, y) {
            circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(4), fileAliases[CELL_META$sector.numeric.index], niceFacing = TRUE)
            circos.lines(c(CELL_META$xlim[1], CELL_META$xlim[2]), c((CELL_META$ylim[1] + CELL_META$ycenter) * 0.5, (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
            for (index in 1:length(x)) {
                circos.lines(c(y[index], y[index]), c(CELL_META$ylim[1], (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
            }
            n <- length(x) - 1
            circos.rect(y[2:(n + 1)] - x[2:(n + 1)], rep(CELL_META$ycenter, n), y[2:(n + 1)], rep(CELL_META$ylim[2], n), col = countRamp(x[2:(n + 1)]), border = NA)
        }
    )

    # Add links
    for (index in 1:nrow(A)) {
        if (A$linkIndexB[index] > 0 && showLinks[1] == TRUE) {
            circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "B", c(B$accumRelCount[A$linkIndexB[index]] - B$relCount[A$linkIndexB[index]], B$accumRelCount[A$linkIndexB[index]]), col = linkRampAB(A$relCount[index] + B$relCount[A$linkIndexB[index]]), border = NA)
        }
        if (A$linkIndexC[index] > 0 && showLinks[2] == TRUE) {
            circos.link("A", c(A$accumRelCount[index] - A$relCount[index], A$accumRelCount[index]), "C", c(C$accumRelCount[A$linkIndexC[index]] - C$relCount[A$linkIndexC[index]], C$accumRelCount[A$linkIndexC[index]]), col = linkRampAC(A$relCount[index] + C$relCount[A$linkIndexC[index]]), border = NA)
        }
    }
    for (index in 1:nrow(B)) {
        if (B$linkIndexC[index] > 0 && showLinks[3] == TRUE) {
            circos.link("B", c(B$accumRelCount[index] - B$relCount[index], B$accumRelCount[index]), "C", c(C$accumRelCount[B$linkIndexC[index]] - C$relCount[B$linkIndexC[index]], C$accumRelCount[B$linkIndexC[index]]), col = linkRampBC(B$relCount[index] + C$relCount[B$linkIndexC[index]]), border = NA)
        }
    }

    # Close diagramm and output file
    circos.clear()
    dev.off()

    # Store parameters for later reference
    parameterFileName <- paste(saveFolder, "circos3.csv", sep = "")
    parameter <- data.frame(
        diagrammName = c(diagrammName),
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
        showLinksBC = c(showLinks[3])
    )
    if (file.exists(parameterFileName) == TRUE) {
        parameters <- read.csv(parameterFileName)
        if (nrow(parameters[parameters$diagrammName == diagrammName, ]) == 0) {
            parameters <- rbind(parameters, parameter)
        } else {
            previousIndex <- match(c(diagrammName), parameters$diagrammName)
            parameters[previousIndex, ] <- parameter
        }
        write.csv(parameters, parameterFileName, row.names = FALSE)
    } else {
        write.csv(parameter, parameterFileName, row.names = FALSE)
    }
}
