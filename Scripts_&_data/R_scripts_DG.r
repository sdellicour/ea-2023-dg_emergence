library(ape)
library(diagram)
library(fields)
library(igraph)
library(lubridate)
library(maptools)
library(MetBrewer)
library(RColorBrewer)
library(raster)
library(rgeos)
library(seraphim)
library(vegan)

# A. Analysis of the genesis

segments = c("ALL_12022025","PA_12022025","PB2_12022025")
	# ALL_12022025 = all concatanated segments but PA and PB2; "ALL_120225_4.fas" = "concatenated" alignment without recombinant sequences
	# (n.b.: the sequence EPI_ISL_19081847 was included two times and the first version was discarded because associated with less metadata)

	# A.1. Analyses based on the alignment including an extended set of background sequences

for (i in 1:length(segments))
	{
		tab1 = read.csv(paste0("Genesis_of_DG/",gsub("2025","25_1",segments[i]),".csv"), sep=";", head=T)
		txt = scan(paste0("Genesis_of_DG/",gsub("2025","25_1",segments[i]),".fas"), what="", sep="\n", quiet=T)
		for (j in 1:length(txt))
			{
				if (grepl(">",txt[j]))
					{
						index = which(tab1[,"SeqID"]==txt[j]); # print(tab1[index,"Isolate_Id"])
						if (length(index) != 1)
							{
								print(c(i,j))
								txt[j] = unlist(strsplit(txt[j],"\\|"))[1]
							}	else	{
								txt[j] = paste0(">",tab1[index,"Isolate_Id"])
							}
					}
			}
		write(txt, paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".fas"))
		tab2 = tab1[,c("Isolate_Id","Collection_Date","Location")]; colnames(tab2) = c("trait","collection_date","location")
		for (j in 1:dim(tab2)[1]) tab2[j,"location"] = unlist(strsplit(tab2[j,"location"]," / "))[2]
		for (j in 1:dim(tab2)[1]) tab2[j,"location"] = gsub("Korea, Republic of","South Korea",tab2[j,"location"])
		for (j in 1:dim(tab2)[1]) tab2[j,"location"] = gsub("Russian Federation","Russia",tab2[j,"location"])
		write.table(tab2, paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".txt"), row.names=F, quote=F, sep="\t")
	}
# system("bash IQTREE_runs_1.sh") # to do rename/sort the IQTREE outputs
samples_DG = read.table("Genesis_of_DG/All_DG_samples.txt", head=F)[,1]
metadata1 = read.table(paste0("Genesis_of_DG/ALL_120225_2.txt"), sep="\t", head=T)
metadata2 = read.table(paste0("Genesis_of_DG/PA_120225_2.txt"), sep="\t", head=T)
metadata3 = read.table(paste0("Genesis_of_DG/PB2_120225_2.txt"), sep="\t", head=T)
countries = unique(rbind(metadata1, metadata2, metadata3)[,"location"]); countries = countries[order(countries)]
countries_colours = colorRampPalette(brewer.pal(11,"Spectral"))(length(countries))
for (i in 1:length(segments)) # ".tre" files have to be preliminary midpoint rooted and re-organised in FigTree
	{
		locations = read.table(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".txt"), sep="\t", head=T)
		tre = read.nexus(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".tre")); rootHeight = max(nodeHeights(tre))
		pdf(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".pdf"), width=6, height=14); par(oma=c(0,0,0,0), mar=c(0,1,0,1))
		plot(tre, x.lim=c(0,rootHeight+(rootHeight/10)), show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.5, col="gray30", edge.color="gray30")
		red = rgb(222,67,39,255,maxColorValue=255)
		for (j in 1:dim(tre$edge)[1])
			{
				if ((!tre$edge[j,2]%in%tre$edge[,1])&&(tre$tip.label[tre$edge[j,2]]%in%samples_DG))
					{
						nodelabels(node=tre$edge[j,2], pch=16, cex=0.40, col=red)
						nodelabels(node=tre$edge[j,2], pch=1, cex=0.40, col="gray30", lwd=0.4)
					}
				if (!tre$edge[j,2]%in%tre$edge[,1])
					{
						country = locations[which(locations[,"trait"]==tre$tip.label[tre$edge[j,2]]),"location"]
						colour = countries_colours[which(countries==country)]
						if (tre$tip.label[tre$edge[j,2]]%in%samples_DG)
							{
								tiplabels(paste0("          ",tre$tip.label[tre$edge[j,2]]), tre$edge[j,2], frame="none", bg=NULL, cex=0.20, col=colour, offset=0.0011)
							}	else	{
								tiplabels(paste0("          ",tre$tip.label[tre$edge[j,2]]), tre$edge[j,2], frame="none", bg=NULL, cex=0.20, col=colour, offset=0.0010)
							}
					}
			}
		if (i == 1)
			{
				add.scale.bar(x=0.038, y=1.1, length=NULL, ask=F, lwd=0.5 , lcol ="gray30", cex=0.6)
				legend(x=0.043, y=580, countries, text.col="gray30", pch=16, pt.cex=1.0, col=countries_colours, box.lty=0, cex=0.5, y.intersp=1.1)
				legend(x=0.043, y=580, countries, text.col=rgb(0,0,0,0), pch=1, pt.cex=1.0, col="gray30", box.lty=0, cex=0.5, pt.lwd=0.3, y.intersp=1.1)
			}
		dev.off()
	}

	# A.1. Analyses based on the alignment only including a restricted set of background sequences

		# A.2.1. Preparation of the metadata files for the discrete phylogeographic analyses

tab2 = read.table(paste0("Genesis_of_DG/ALL_120225_2.txt"), head=T, sep="\t")
selected_clade = read.table(paste0("Genesis_of_DG/Selected_clade.txt"), head=F, sep="\t")
	# see the clade highlighted in "ALL_120225_2.png" by a red arroiw
tab3 = tab2[which(tab2[,"trait"]%in%selected_clade[,1]),]
write.table(tab3, "Genesis_of_DG/ALL_120225_3.txt", row.names=F, quote=F, sep="\t")
tab3 = read.table(paste0("Genesis_of_DG/ALL_120225_3.txt"), head=T, sep="\t"); txt3 = c()
txt2 = scan(paste0("Genesis_of_DG/ALL_120225_2.fas"), what="", sep="\n", quiet=T); tmp1 = c()
for (i in 1:length(txt2))
	{
		if (grepl(">",txt2[i]))
			{
				tmp1 = c(tmp1, txt2[i])
			}	else	{
				if (grepl(">",txt2[i-1]))
					{
						tmp2 = txt2[i]
					}
				if (!grepl(">",txt2[i-1]))
					{
						tmp2 = paste0(tmp2,txt2[i])
					}
				if ((grepl(">",txt2[i+1]))|(i == length(txt2)))
					{
						tmp1 = c(tmp1, tmp2)
						tmp2 = NULL
					}
			}
	}
txt2 = tmp1
for (i in 1:length(txt2))
	{
		if (grepl(">",txt2[i]))
			{
				if (gsub(">","",txt2[i])%in%tab3[,"trait"])
					{
						txt3 = c(txt3, txt2[i], txt2[i+1])
					}
			}
	}
write(txt3, paste0("Genesis_of_DG/ALL_120225_3.fas"))
	
		# A.2.2. Conducting the recombination analyses with the phi-test and RDP4

txt4 = scan(paste0("Genesis_of_DG/ALL_120225_4.fas"), what="", sep="\n", quiet=T)
	# corresponds to All_RDP_analyses/ALL_120225_3_res/ALL_120225_3_without_recombinant_sequences.fas
tab3 = read.table(paste0("Genesis_of_DG/ALL_120225_3.txt"), head=T, sep="\t")
selected_sequences = gsub(">","",txt4[which(grepl(">",txt4))])
samples_DG = read.table("Genesis_of_DG/All_DG_samples.txt", head=F)[,1]
print(samples_DG[which(!samples_DG%in%selected_sequences)]) # none - ok
tab4 = tab3[which(tab3[,"trait"]%in%selected_sequences),]
write.table(tab4, "Genesis_of_DG/ALL_120225_4.txt", row.names=F, quote=F, sep="\t")

		# A.2.3. Retrieving and annotating the MCC tree for each segment BEAST analysis

burnIns = rep(NA, length(segments))
for (i in 1:length(segments))
	{
		if (segments[i] == "ALL_12022025") fileName = paste0("Genesis_of_DG/",gsub("2025","25_4",segments[i]),".log")
		if (segments[i] != "ALL_12022025") fileName = paste0("Genesis_of_DG/",gsub("2025","25_2",segments[i]),".log")		
		txt = scan(fileName, what="", sep="\n", quiet=T); txt = txt[length(txt)]
		nberOfSampledTrees = as.numeric(unlist(strsplit(txt,"\t"))[1])/100000; burnIns[i] = (nberOfSampledTrees/10)+1
	}
for (i in 1:length(segments))
	{
		if (segments[i] == "ALL_12022025")
			{
				system(paste0("Genesis_of_DG/BEAST_v_1_10_4/bin/treeannotator -burninTrees ",round(burnIns[i])," -heights keep Genesis_of_DG/",gsub("2025","25_4",segments[i]),".trees Genesis_of_DG/",gsub("2025","25_4",segments[i]),"_NEW.tree"), ignore.stdout=F, ignore.stderr=F) # MCC trees to be re-organised with CMD+U in FigTree
			}	else	{
				system(paste0("Genesis_of_DG/BEAST_v_1_10_4/bin/treeannotator -burninTrees ",round(burnIns[i])," -heights keep Genesis_of_DG/",gsub("2025","25_2",segments[i]),".trees Genesis_of_DG/",gsub("2025","25_2",segments[i]),"_NEW.tree"), ignore.stdout=F, ignore.stderr=F) # MCC trees to be re-organised with CMD+U in FigTree				
			}
	}
		# Information retrieved from the three MCC trees (TO BE UPDATED):
			# - most probable ancestral location inferred for the most recent common ancestor of the DG samples clade: Sweden (concatenated minus PA and PB2,
			#   with a posterior probability = 0.998), and Germany (for both PA and PB2, with a posterior probability >0.999 for both)
			# - time of the most recent common ancestor of the DG samples clade: 2023.65 (concatenated minus PA and PB2, 95% HPD = [2023.56-2023.73]), 
			#   2023.60 (PA, 95% HPD = [2023.41-2023.77]), and 2023.12 (PB2, 95% HPD = [2022.90-2023.38])

		# A.2.3. Visualisation of the discrete phylogeographic reconstruction per segment

samples_DG = read.table("Genesis_of_DG/All_DG_samples.txt", head=F)[,1]
for (h in 1:length(segments))
	{
		if (segments[h] == "ALL_12022025")
			{
				tab3 = read.table(paste0("Genesis_of_DG/",gsub("2025","25_4",segments[h]),".txt"), head=T, sep="\t")
				mcc_tre = readAnnotatedNexus(paste0("Genesis_of_DG/",gsub("2025","25_4",segments[h]),".tree"))
			}	else	{
				tab3 = read.table(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[h]),".txt"), head=T, sep="\t")
				mcc_tre = readAnnotatedNexus(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[h]),".tree"))				
			}
		for (i in 0:length(mcc_tre$annotations))
			{
				if (i == 0)
					{
						if (mcc_tre$root.annotation$location == "Korea") mcc_tre$root.annotation$location = "North Korea"
						if (mcc_tre$root.annotation$location == "DemocraticPeople'sRepublicof") mcc_tre$root.annotation$location = "North Korea"
						if (mcc_tre$root.annotation$location == "Korea,DemocraticPeople'sRepublicof") mcc_tre$root.annotation$location = "North Korea"
						for (j in 1:length(mcc_tre$root.annotation$location.set))
							{
								if (mcc_tre$root.annotation$location.set[[j]] == "Korea") mcc_tre$root.annotation$location.set[[j]] = "North Korea"
								if (mcc_tre$root.annotation$location.set[[j]] == "DemocraticPeople'sRepublicof") mcc_tre$root.annotation$location.set[[j]] = "North Korea"
							}						
					}
				if (i != 0)
					{
						if (mcc_tre$annotations[[i]]$location == "Korea") mcc_tre$annotations[[i]]$location = "North Korea"
						if (mcc_tre$annotations[[i]]$location == "DemocraticPeople'sRepublicof") mcc_tre$annotations[[i]]$location = "North Korea"
						if (mcc_tre$annotations[[i]]$location == "Korea,DemocraticPeople'sRepublicof") mcc_tre$annotations[[i]]$location = "North Korea"
						for (j in 1:length(mcc_tre$annotations[[i]]$location.set))
							{
								if (mcc_tre$annotations[[i]]$location.set[[j]] == "Korea") mcc_tre$annotations[[i]]$location.set[[j]] = "North Korea"
								if (mcc_tre$annotations[[i]]$location.set[[j]] == "DemocraticPeople'sRepublicof") mcc_tre$annotations[[i]]$location.set[[j]] = "North Korea"
							}
					}
			}
		collection_dates = decimal_date(ymd(tab3[,"collection_date"]))
		mostRecentSamplingDatum = max(collection_dates, na.rm=T)
		rootHeight = max(nodeHeights(mcc_tre)); root_time = mostRecentSamplingDatum-rootHeight
		minYear = mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]]
		maxYear = mostRecentSamplingDatum; nodes_cols = c(); countries = c()
		for (i in 1:length(mcc_tre$annotations)) countries = c(countries, mcc_tre$annotations[[i]]$location)
		countries = unique(countries); countries = countries[order(countries)]
		colours = colorRampPalette(brewer.pal(11,"Spectral"))(length(countries))
		root_countries_prob = matrix(0, nrow=1, ncol=length(countries))
		countries_prob = matrix(0, nrow=length(mcc_tre$annotations), ncol=length(countries))
		colnames(root_countries_prob) = countries; colnames(countries_prob) = countries
		index = which(countries == mcc_tre$root.annotation$location); root_node_col = colours[index]
		for (i in 1:length(mcc_tre$root.annotation$location.set))
			{
				root_countries_prob[1,mcc_tre$root.annotation$location.set[[i]][1]] = mcc_tre$root.annotation$location.set.prob[[i]][1]
			}
		for (i in 1:length(mcc_tre$annotations))
			{
				index = which(countries == mcc_tre$annotations[[i]]$location); nodes_cols = c(nodes_cols, colours[index])
				for (j in 1:length(mcc_tre$annotations[[i]]$location.set))
					{
						countries_prob[i,mcc_tre$annotations[[i]]$location.set[[j]][1]] = mcc_tre$annotations[[i]]$location.set.prob[[j]][1]
					}
			}
		pieCharts = TRUE; tipNodes = c()
		for (i in 1:dim(mcc_tre$edge)[1])
			{
				if (!mcc_tre$edge[i,2]%in%mcc_tre$edge[,1]) tipNodes = c(tipNodes, TRUE)
				if (mcc_tre$edge[i,2]%in%mcc_tre$edge[,1]) tipNodes = c(tipNodes, FALSE)
			}
		root_countries_prob[root_countries_prob[]<0.05] = 0; countries_prob[countries_prob[]<0.05] = 0
		if (segments[h] == "ALL_12022025")
			{
				pdf(paste0("Genesis_of_DG/",gsub("2025","25_4",segments[h]),"_NEW.pdf"), width=8, height=4); # dev.new(width=8, height=4)
			}	else	{
				pdf(paste0("Genesis_of_DG/",gsub("2025","25_2",segments[h]),"_NEW.pdf"), width=8, height=4); # dev.new(width=8, height=4)
			}
		par(mar=c(0.2,2,0,0), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = TRUE
		plot(mcc_tre, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, direction="downwards",
			 y.lim=c(0, max(nodeHeights(mcc_tre))), col="gray30", edge.color="gray30")
		mcc_tre_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); rootBarPlotted = FALSE
		gray90_transparent = rgb(229, 229, 229, 150, names=NULL, maxColorValue=255)
		for (j in 1:dim(mcc_tre$edge)[1])
			{
				if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&(length(mcc_tre$annotations[[j]]$`height_95%_HPD`) > 1))
					{
						y1 = mostRecentSamplingDatum-(mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[2]])
						y2 = mostRecentSamplingDatum-(mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[1]])
						lines(x=rep(mcc_tre_obj$xx[mcc_tre$edge[j,2]],2), y=c(y2,y1), lwd=3.5, lend=0, col=gray90_transparent)
					}
				if ((rootBarPlotted == FALSE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
					{
						y1 = mostRecentSamplingDatum-(mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]])
						y2 = mostRecentSamplingDatum-(mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[1]])
						lines(x=rep(mcc_tre_obj$xx[mcc_tre$edge[j,1]],2), y=c(y2,y1), lwd=3.5, lend=0, col=gray90_transparent)
						rootBarPlotted = TRUE
					}
			}
		plotted_root = FALSE
		for (i in 1:dim(mcc_tre$edge)[1])
			{
				if (pieCharts == FALSE)
					{
						nodelabels(node=mcc_tre$edge[i,2], pch=16, cex=0.8, col=nodes_cols[i])
						nodelabels(node=mcc_tre$edge[i,2], pch=1, cex=0.8, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
					}	else	{
						if ((sum(countries_prob[i,]==1) > 0)|(sum(countries_prob[i,]==0) == (length(countries)-1)))
							{
								nodelabels(node=mcc_tre$edge[i,2], cex=0.15, pie=t(countries_prob[i,]), piecol=colours)
							}	else		{
								nodelabels(node=mcc_tre$edge[i,2], cex=0.45, pie=t(countries_prob[i,]), piecol=colours)
							}	
					}
				if ((plotted_root == FALSE)&&(!mcc_tre$edge[i,1]%in%mcc_tre$edge[,2]))
					{
						if (pieCharts == FALSE)
							{
								nodelabels(node=mcc_tre$edge[i,1], pch=16, cex=0.8, col=nodes_cols[i])
								nodelabels(node=mcc_tre$edge[i,1], pch=1, cex=0.8, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
							}	else	{
								if ((sum(root_countries_prob[i,]==1) > 0)|(sum(root_countries_prob[i,]==0) == (length(countries)-1)))
									{
										nodelabels(node=mcc_tre$edge[i,1], cex=0.15, pie=t(root_countries_prob[i,]), piecol=colours)
									}	else		{
										nodelabels(node=mcc_tre$edge[i,1], cex=0.45, pie=t(root_countries_prob[i,]), piecol=colours)
									}	
							}
						plotted_root = TRUE
					}
			}
		for (i in 1:dim(mcc_tre$edge)[1])
			{
				if (!mcc_tre$edge[i,2]%in%mcc_tre$edge[,1])
					{
						if (mcc_tre$tip.label[mcc_tre$edge[i,2]]%in%samples_DG)
							{
								# tiplabels(paste0(mcc_tre$tip.label[mcc_tre$edge[i,2]]), mcc_tre$edge[i,2], frame="none", bg=NULL, cex=0.20, col="gray35", offset=0.0011)
								tiplabels("X", mcc_tre$edge[i,2], frame="none", bg=NULL, cex=0.30, col="black", offset=0.05, asp=c(0,2))
							}
					}
			}
		if (h == 1)
			{
				abline(h=mostRecentSamplingDatum-2022, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2023, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2024, lwd=1, col="gray90", lty=2)
				dates_to_plot = c(2019,2020,2021,2022,2023,2024,2025); dates_to_print = dates_to_plot
				axis(lwd=0.3, at=mostRecentSamplingDatum-dates_to_plot, labels=dates_to_print, cex.axis=0.60, 
					 mgp=c(0,0.15,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.013, side=2, las=1)
				legend(x=155, y=2.5, countries, text.col="gray30", pch=16, pt.cex=1.2, col=colours, box.lty=0, cex=0.6, x.intersp=0.8, y.intersp=1.1)
				legend(x=155, y=2.5, countries, text.col=rgb(0,0,0,0), pch=1, pt.cex=1.2, col="gray30", box.lty=0, cex=0.6, pt.lwd=0.3, x.intersp=0.8, y.intersp=1.1)
			}
		if (h == 2)
			{
				abline(h=mostRecentSamplingDatum-2010, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2015, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2020, lwd=1, col="gray90", lty=2)
				dates_to_plot = seq(2005,2035,5); dates_to_print = dates_to_plot
				axis(lwd=0.3, at=mostRecentSamplingDatum-dates_to_plot, labels=dates_to_print, cex.axis=0.60, 
					 mgp=c(0,0.15,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.013, side=2, las=1)
				legend(x=450, y=17, countries, text.col="gray30", pch=16, pt.cex=1.2, col=colours, box.lty=0, cex=0.6, x.intersp=0.8, y.intersp=1.1)
				legend(x=450, y=17, countries, text.col=rgb(0,0,0,0), pch=1, pt.cex=1.2, col="gray30", box.lty=0, cex=0.6, pt.lwd=0.3, x.intersp=0.8, y.intersp=1.1)
			}
		if (h == 3)
			{
				abline(h=mostRecentSamplingDatum-1990, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2000, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2010, lwd=1, col="gray90", lty=2)
				abline(h=mostRecentSamplingDatum-2020, lwd=1, col="gray90", lty=2)
				dates_to_plot = seq(1980,2030,10); dates_to_print = dates_to_plot
				axis(lwd=0.3, at=mostRecentSamplingDatum-dates_to_plot, labels=dates_to_print, cex.axis=0.60, 
					 mgp=c(0,0.15,-0.3), lwd.tick=0.5, col.lab="gray30", col="gray30", tck=-0.013, side=2, las=1)
				legend(x=455, y=35, countries, text.col="gray30", pch=16, pt.cex=1.2, col=colours, box.lty=0, cex=0.6, x.intersp=0.8, y.intersp=1.1)
				legend(x=455, y=35, countries, text.col=rgb(0,0,0,0), pch=1, pt.cex=1.2, col="gray30", box.lty=0, cex=0.6, pt.lwd=0.3, x.intersp=0.8, y.intersp=1.1)
			}			
		dev.off()
	}

# B. Analysis of the spread

localTreesDirectory = "Spread_of_DG/ALL_02122024_ext"; nberOfExtractionFiles = 1000
collection_dates = read.table("Spread_of_DG/ALL_021224_2.txt", head=T)[,"collection_date"]
mostRecentSamplingDatum = max(decimal_date(ymd(collection_dates)))

	# B.1. Extracting the spatio-temporal information embedded in posterior trees

source("Extractions_1.r") # for the MCC tree
source("Extractions_2.r") # for the posterior trees
nberOfTreesToSample = nberOfExtractionFiles; burnIn = 251
allTrees1 = readAnnotatedNexus("Spread_of_DG/ALL_021224_2.trees")
allTrees2 = allTrees1[(burnIn+1):length(allTrees1)]
for (i in 1:length(allTrees2))
	{
		csv = Extractions_2(allTrees2[[i]], mostRecentSamplingDatum)
		write.csv(csv, paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), row.names=F, quote=F)
	}
mcc_tre = readAnnotatedNexus("Spread_of_DG/ALL_021224_2.tree")
mcc_tab = Extractions_1(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "Spread_of_DG/ALL_021224_2.csv", row.names=F, quote=F)

	# B.2. Estimating the dispersal statistics associated with H5N1 DG lineages

timeSlices = 100; onlyTipBranches = F; showingPlots = F; nberOfCores = 10; slidingWindow = 1/24
outputName = paste0("Spread_of_DG/Dispersal_statistics/ALL_021224")
spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)

mat = read.table("Spread_of_DG/Dispersal_statistics/ALL_021224_estimated_dispersal_statistics.txt", head=T)
vS1 = mat[,"weighted_diffusion_coefficient"]; HPD1 = round(HDInterval::hdi(vS1)[1:2],0)
vS2 = mat[,"isolation_by_distance_signal_rP2"]; HPD2 = round(HDInterval::hdi(vS2)[1:2],3)
cat("WDC = ",round(median(vS1),0)," km2/year (95% HPD = [",HPD1[1],", ",HPD1[2],"])",sep="")
	# WDC = XXXX km2/year (95% HPD = [XXXX, XXXX])
cat("IBD (rP2) = ",round(median(vS2),3)," (95% HPD = [",HPD2[1],", ",HPD2[2],"])",sep="")
	# IBD (rP2) = XXXX (95% HPD = [XXXX, XXXX])

	# B.3. Mapping the inferred dispersal history of inferred H5N1 DG lineages

e_Palearctic = extent(-13, 150, 35, 73); e_Europe = extent(-13, 55, 35, 73)
e_studyArea_1 = extent(-11, 29, 41.5, 65); e_studyArea_2 = extent(-3, 22, 45, 61); e_studyArea_3 = extent(-3, 22, 46, 62)
countries1 = crop(shapefile("Spread_of_DG/Different_shapefiles/World_countries_shps/World_countries_shapefile.shp"), e_Palearctic)
borders1 = crop(shapefile("Spread_of_DG/Different_shapefiles/International_borders/Only_international_borders.shp"), e_Palearctic)
coasts1 = crop(shapefile("Spread_of_DG/Different_shapefiles/Coast_lines_borders/Only_coast_lines_borders.shp"), e_Palearctic)
countries2 = crop(countries1, e_Europe); borders2 = crop(borders1, e_Europe); coasts2 = crop(coasts1, e_Europe)
countries3 = crop(countries1, e_studyArea_1); borders3 = crop(borders1, e_studyArea_1); coasts3 = crop(coasts1, e_studyArea_1)
countries4 = crop(countries1, e_studyArea_2); borders4 = crop(borders1, e_studyArea_2); coasts4 = crop(coasts1, e_studyArea_2)
countries5 = crop(countries1, e_studyArea_3); borders5 = crop(borders1, e_studyArea_3); coasts5 = crop(coasts1, e_studyArea_3)

mcc = read.csv("Spread_of_DG/ALL_021224_2.csv", head=T)
mcc = mcc[order(mcc[,"startYear"]),]; mcc1 = mcc[1,]; mcc2 = mcc[c(2:dim(mcc)[1]),]
mcc2 = mcc2[order(mcc2[,"endYear"]),]; mcc = rbind(mcc1,mcc2)
mcc_tre = readAnnotatedNexus("Spread_of_DG/ALL_021224_2.tree"); mcc_tre$tip.label = gsub("'","",mcc_tre$tip.label)
rootHeight = max(nodeHeights(mcc_tre)); root_time = mostRecentSamplingDatum-rootHeight
minYear = mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]]; maxYear = mostRecentSamplingDatum
hpd95_rootHeight = c(mcc_tre$root.annotation$`height_95%_HPD`[[2]],mcc_tre$root.annotation$`height_95%_HPD`[[1]])
print(c(root_time,mostRecentSamplingDatum-hpd95_rootHeight)) # root time = 2023.59 (95% HPD = [2023.46, 2023.68])
selectedLabels = c("2023-07-01","2023-09-01","2023-11-01","2024-01-01","2024-03-01","2024-05-01")
selectedDates = c(minYear, decimal_date(ymd(selectedLabels)), maxYear); selectedLabels = c("", selectedLabels, "")
colour_scale = met.brewer(name="Hiroshige", n=111, type="continuous")[1:101]

prob = 0.95; precision = 1/24; startDatum = minYear
polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))
endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
	{
		date = as.numeric(names(polygons[[i]]))
		polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
		polygons_colours[i] = paste0(colour_scale[polygon_index],"20")
	}

pdf("Spread_of_DG/ALL_021224_2_NEW1.pdf", width=4.0, height=4.28); # dev.new(width=4.0, height=4.05)
par(mar=c(0.8,0,0,0.01), oma=c(0,0,0,0), mgp=c(0,0.1,0), lwd=0.2, bty="o", col="gray30"); plottingRootNode = TRUE
plot(mcc_tre, show.tip.label=F, show.node.label=F, edge.width=0.5, cex=0.6, align.tip.label=3, 
	 x.lim=c(minYear-(maxYear-max(nodeHeights(mcc_tre))), max(nodeHeights(mcc_tre))+0.085), col="gray30", edge.color="gray30")
mcc_tre_obj = get("last_plot.phylo", envir=.PlotPhyloEnv); rootBarPlotted = FALSE
for (j in 1:dim(mcc_tre$edge)[1])
	{
		endYear = root_time+nodeHeights(mcc_tre)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&(length(mcc_tre$annotations[[j]]$`height_95%_HPD`) > 1))
			{
				x1 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$annotations[[j]]$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,2]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
			}
		if ((rootBarPlotted == FALSE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
			{
				endYear = root_time+nodeHeights(mcc_tre)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]
				x1 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[2]])-root_time
				x2 = (mostRecentSamplingDatum-mcc_tre$root.annotation$`height_95%_HPD`[[1]])-root_time
				lines(x=c(x1,x2), y=rep(mcc_tre_obj$yy[mcc_tre$edge[j,1]],2), lwd=3.5, lend=0, col=paste0(endYear_colour,"40"))
				rootBarPlotted = TRUE
			}				
	}
for (j in 1:dim(mcc_tre$edge)[1])
	{
		endYear = root_time+nodeHeights(mcc_tre)[j,2]
		endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
		endYear_colour = colour_scale[endYear_index]
		if ((mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])&&(mcc_tre$annotations[[j]]$posterior >= 0.95))
			{
				nodelabels(node=mcc_tre$edge[j,2], pch=16, cex=0.70, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,2], pch=1, cex=0.70, col="gray30", lwd=0.2)
			}
		if (!mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])
			{
				nodelabels(node=mcc_tre$edge[j,2], pch=15, cex=0.55, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,2], pch=0, cex=0.55, col="gray30", lwd=0.2)
			}
		if ((plottingRootNode == TRUE)&&(!mcc_tre$edge[j,1]%in%mcc_tre$edge[,2]))
			{
				endYear = root_time+nodeHeights(mcc_tre)[j,1]
				endYear_index = (((endYear-minYear)/(maxYear-minYear))*100)+1
				endYear_colour = colour_scale[endYear_index]; plottingRootNode = TRUE		
				nodelabels(node=mcc_tre$edge[j,1], pch=16, cex=0.70, col=endYear_colour)
				nodelabels(node=mcc_tre$edge[j,1], pch=1, cex=0.70, col="gray30", lwd=0.2)
			}
	}
for (j in 1:dim(mcc_tre$edge)[1])
	{
		if (!mcc_tre$edge[j,2]%in%mcc_tre$edge[,1])
			{
				tipLabel = paste0("   ",mcc_tre$tip.label[mcc_tre$edge[j,2]])
				tiplabels(tipLabel, tip=mcc_tre$edge[j,2], col="gray30", frame="none", bg=NULL, cex=0.35, adj=c(0.0,0.5))
			}
	}
axis(lwd=0.3, at=selectedDates-root_time, labels=selectedLabels, cex=0.40, cex.axis=0.40, mgp=c(0,-0.3,-0.3),
	 lwd.tick=0.3, col.axis="gray30", col.lab="gray30", col="gray30", tck=-0.008, side=1)
dev.off()

selectedLabels = c("2023-09-01","2023-11-01","2024-01-01"); selectedDates = decimal_date(ymd(selectedLabels))
cutOffs = c(selectedDates, maxYear); croppingPolygons = FALSE; plottingAllNodes = FALSE
for (h in length(cutOffs))
	{
		pdf(paste0("Spread_of_DG/ALL_021224_2_NEW2.pdf"), width=5.2, height=5)
		par(oma=c(0,0,0,0), mar=c(0.8,0.7,1,0.4), mgp=c(0,0.1,0), lwd=0.2, bty="o")
		plot(countries3, col="gray90", border=NA, ann=F, axes=F)
		plot(borders3, col="white", lwd=0.3, add=T)
		plot(coasts3, col="gray70", lwd=0.5, add=T)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		if (h == length(cutOffs))
			{
				plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.10,0.40,0.895,0.905),
	 				 legend.args=list(text="", cex=0.7, col="gray30"), horizontal=T,
					 axis.args=list(cex.axis=0.45, lwd=0, lwd.tick=0.2, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0, mgp=c(0,-0.20,0)))
			 }
		if (plottingAllNodes == TRUE)
			{
				for (i in 1:nberOfExtractionFiles)
					{
						csv = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
						internalNodes = which(csv[,"node2"]%in%csv[,"node1"]); tipNodes = which(!csv[,"node2"]%in%csv[,"node1"])
						points(csv[internalNodes,c("endLon","endLat")], cex=0.5, pch=1, lwd=0.1, col="gray30")
						if (i == nberOfExtractionFiles)
							{
								points(csv[tipNodes,c("endLon","endLat")], cex=0.3, pch=16, col="red")
							}
					}
				internalNodes = which(mcc[,"node2"]%in%mcc[,"node1"])
				points(mcc[internalNodes,c("endLon","endLat")], cex=0.3, pch=16, col="black")
			}
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								# polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(countries2)
						if (croppingPolygons == TRUE) pol = crop(pol, countries2)
						plot(pol, axes=F, col=polygons_colours[i], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in dim(mcc)[1]:1)
			{
				if ((mcc[i,"startYear"] <= cutOffs[h])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.6)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.3, cex=0.6)
					}
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						if (mcc[i,"node2"]%in%mcc[,"node1"])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.6)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.3, cex=0.6)
							}	else	{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.5)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.3, cex=0.5)					
							}
					}
			}
		rect(-11, 41.5, 29, 65, lwd=0.2, border="gray30")
		dev.off()
	}

selectedLabels = c("2023-09-01","2023-11-01","2024-01-01"); selectedDates = decimal_date(ymd(selectedLabels))
cutOffs = c(selectedDates, maxYear); croppingPolygons = FALSE; plottingAllNodes = FALSE
pdf(paste0("Spread_of_DG/ALL_021224_2_NEW3.pdf"), width=11.5, height=3)
par(mfrow=c(1,4), oma=c(0,0,0.1,0), mar=c(0.1,0.1,0.1,0.1), mgp=c(0,0.1,0), lwd=0.2, bty="o")
selectedLabels = c(selectedLabels, as.character(max(ymd(collection_dates))))
for (h in 1:length(cutOffs))
	{
		plot(countries4, col="gray90", border=NA, ann=F, axes=F)
		plot(borders4, col="white", lwd=0.3, add=T)
		plot(coasts4, col="gray70", lwd=0.5, add=T)
		rast = raster(matrix(nrow=1, ncol=2)); rast[1] = startDatum; rast[2] = max(mcc[,"endYear"])
		if (plottingAllNodes == TRUE)
			{
				for (i in 1:nberOfExtractionFiles)
					{
						csv = read.csv(paste0(localTreesDirectory,"/TreeExtractions_",i,".csv"), head=T)
						internalNodes = which(csv[,"node2"]%in%csv[,"node1"]); tipNodes = which(!csv[,"node2"]%in%csv[,"node1"])
						points(csv[internalNodes,c("endLon","endLat")], cex=0.5, pch=1, lwd=0.1, col="gray30")
						if (i == nberOfExtractionFiles)
							{
								points(csv[tipNodes,c("endLon","endLat")], cex=0.3, pch=16, col="red")
							}
					}
				internalNodes = which(mcc[,"node2"]%in%mcc[,"node1"])
				points(mcc[internalNodes,c("endLon","endLat")], cex=0.3, pch=16, col="black")
			}
		for (i in 1:length(polygons))
			{
				if (as.numeric(names(polygons[[i]])) <= cutOffs[h])
					{
						for (j in 1:length(polygons[[i]]@polygons))
							{
								# polygons[[i]]@polygons[[j]] = checkPolygonsHoles(polygons[[i]]@polygons[[j]])
							}
						pol = polygons[[i]]; crs(pol) = crs(countries2)
						if (croppingPolygons == TRUE) pol = crop(pol, countries2)
						plot(pol, axes=F, col=polygons_colours[[i]][j], add=T, border=NA)
					}
			}
		for (i in 1:dim(mcc)[1])
			{
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
				  		  		    arr.width=0, lwd=0.3, lty=1, lcol="gray30", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
				  	}
			}
		for (i in dim(mcc)[1]:1)
			{
				if ((mcc[i,"startYear"] <= cutOffs[h])&(!mcc[i,"node1"]%in%mcc[,"node2"]))
					{
						startYears_index = (((mcc[i,"startYear"]-minYear)/(maxYear-minYear))*100)+1
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=colour_scale[startYears_index], cex=0.9)
						points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray30", lwd=0.3, cex=0.9)
					}
				if (mcc[i,"endYear"] <= cutOffs[h])
					{
						if (mcc[i,"node2"]%in%mcc[,"node1"])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=0.9)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray30", lwd=0.3, cex=0.9)
							}	else	{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=15, col=endYears_colours[i], cex=0.8)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=0, col="gray30", lwd=0.3, cex=0.8)					
							}
					}
			}
		rect(-2.9, 59.5, 4, 60.9, lwd=0.0, border=NA, col="white")
		# rect(-2.9, 60.5, 4, 61.9, lwd=0.0, border=NA, col="white")
		mtext(selectedLabels[h], cex=0.6, col="gray30", at=0.8, line=-2.5)
		rect(-3, 45, 22, 61, lwd=0.2, border="gray30")
		# rect(-3, 46, 22, 62, lwd=0.2, border="gray30")
	}
dev.off()

