
library(xlsx)
library(circlize)

cols = c("dodgerblue2","#E31A1C", # red
         "green4",
         "#6A3D9A", # purple
         "#FF7F00", # orange
         "gold1",
         "skyblue2","#FB9A99", # lt pink
         "palegreen2",
         "#CAB2D6", # lt purple
         "#FDBF6F", # lt orange
         "khaki2",
         "maroon","orchid1","deeppink1","blue1","steelblue4",
         "darkturquoise","green1","yellow4","yellow3",
         "darkorange4","brown")

# must be run inside an analysis directory

if (!dir.exists("chords"))
{
    dir.create("chords", mode="0755", recursive=T)
}

if (!dir.exists("barplots"))
{
    dir.create("barplots", mode="0755", recursive=T)
}

trav_names = c()
df_list = list()

for (f in Sys.glob("full/mid*"))
{
    name = gsub(".csv$", "", basename(f))

    d=read.csv(f, sep="\t")
    d$bestVFamily = gsub("([DN]|/DV12)$", "", d$bestVFamily)

    families = d$bestVFamily
    table_f = as.data.frame(table(families))
    table_f$num = as.numeric(gsub("TR[AB]V","",table_f$families))
    table_f = table_f[order(table_f$num), ]
    print(table_f)

    svg(paste0("barplots/", name, ".svg"))
    par(mar=c(4,8,2,2))
    barplot(table_f$Freq, names.arg=table_f$families, horiz=T, las=2, cex.names=0.7, col=cols[table_f$num])
    dev.off()

    df = data.frame(aggregate(cloneFraction ~ bestVFamily, data=d, sum), row.names=1)
    df$num = as.numeric(gsub("TR[AB]V","", rownames(df)))
    df = df[order(df$num),]
    print(df)

    svg(paste0("barplots/", name, "_abund.svg"))
    par(mar=c(4,8,2,2))
    barplot(df$cloneFraction, names.arg=rownames(df), horiz=T, las=2, cex.names=0.7, col=cols[df$num])
    dev.off()

    trav_names = c(trav_names, rownames(df))

    df_list[[toupper(gsub("_clones", "", name))]] = df

    svg(paste0("chords/", name, "_chord.svg"))
    d$bestVFamily = gsub("TR[AB]", "", d$bestVFamily)
    d$bestJFamily = gsub("TR[AB]", "", d$bestJFamily)
    tab = table(d$bestVFamily, d$bestJFamily)
    chorddf = as.data.frame(tab)
    grid.cols = c(rep("grey",length(unique(chorddf$Var2))), rainbow(length(unique(chorddf$Var1))))
    labels = c(names(sort(colSums(tab), decreasing=T)), names(sort(rowSums(tab), decreasing=T)))
    chordDiagram(chorddf, link.rank=rank(chorddf[[3]]), order=labels, transparency=0.05, grid.col=grid.cols, preAllocateTracks=list(track.height=max(strwidth(unlist(dimnames(chorddf))))), link.sort=T, link.decreasing=F, annotationTrack="grid", annotationTrackHeight=0.03)
    circos.track(track.index=1, panel.fun=function(x, y) {circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, cex=0.5, facing="clockwise", niceFacing=T, adj=c(0, 0.5))}, bg.border=NA)
    dev.off()
}

trav_names = unique(trav_names)
nums = as.numeric(gsub("TR[AB]V","", trav_names))
trav_names = trav_names[order(nums)]
df = data.frame(row.names=trav_names)

for (name in names(df_list))
{
    temp=df_list[[name]]
    
    df[rownames(temp), name] = temp$cloneFraction
}

df[is.na(df)] = 0
print(df)

write.xlsx(t(df), file="table.xlsx")

