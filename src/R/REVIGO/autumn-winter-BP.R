

# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0010035","response to inorganic substance",3.699,4.9251,0.889,0.000,"response to inorganic substance"),
c("GO:0009628","response to abiotic stimulus",7.946,0.4983,0.900,0.433,"response to inorganic substance"),
c("GO:0007267","cell-cell signaling",0.371,0.3142,0.861,0.346,"response to inorganic substance"),
c("GO:0001101","response to acid chemical",5.011,3.4933,0.887,0.580,"response to inorganic substance"),
c("GO:1901700","response to oxygen-containing compound",6.504,2.7624,0.885,0.638,"response to inorganic substance"),
c("GO:0009637","response to blue light",0.341,0.4779,0.920,0.223,"response to inorganic substance"),
c("GO:0030522","intracellular receptor signaling pathway",0.078,0.6435,0.826,0.191,"response to inorganic substance"),
c("GO:0042221","response to chemical",12.434,0.6435,0.897,0.377,"response to inorganic substance"),
c("GO:0071495","cellular response to endogenous stimulus",3.988,0.3357,0.905,0.307,"response to inorganic substance"),
c("GO:0010876","lipid localization",0.708,5.7295,0.862,0.000,"lipid localization"),
c("GO:0055085","transmembrane transport",5.175,0.6851,0.863,0.478,"lipid localization"),
c("GO:0070972","protein localization to endoplasmic reticulum",0.194,0.6851,0.863,0.627,"lipid localization"),
c("GO:0008643","carbohydrate transport",0.630,0.3357,0.829,0.554,"lipid localization"),
c("GO:0072657","protein localization to membrane",0.259,0.9255,0.741,0.542,"lipid localization"),
c("GO:0006826","iron ion transport",0.121,0.4576,0.891,0.314,"lipid localization"),
c("GO:0006857","oligopeptide transport",0.155,3.4917,0.881,0.321,"lipid localization"),
c("GO:0015849","organic acid transport",0.682,0.6851,0.829,0.476,"lipid localization"),
c("GO:0051179","localization",12.352,0.4576,0.972,0.000,"localization"),
c("GO:2000060","positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolic process",0.004,300.0000,0.826,0.000,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0009894","regulation of catabolic process",0.414,0.6851,0.853,0.297,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0018195","peptidyl-arginine modification",0.035,0.6435,0.863,0.200,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0051336","regulation of hydrolase activity",0.781,0.6851,0.912,0.120,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0080090","regulation of primary metabolic process",13.681,0.9845,0.809,0.435,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0018212","peptidyl-tyrosine modification",0.121,0.4576,0.854,0.536,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0097264","self proteolysis",0.004,0.3695,0.894,0.272,"positive regulation of protein ubiquitination involved in ubiquitin-dependent protein catabolism"),
c("GO:0045229","external encapsulating structure organization",2.577,5.5273,0.890,0.032,"external encapsulating structure organization"),
c("GO:0007010","cytoskeleton organization",1.032,0.9381,0.896,0.426,"external encapsulating structure organization"),
c("GO:0031163","metallo-sulfur cluster assembly",0.142,0.6464,0.907,0.343,"external encapsulating structure organization"),
c("GO:0008037","cell recognition",0.220,1.0443,0.873,0.045,"cell recognition"),
c("GO:0044042","glucan metabolic process",0.962,0.3597,0.867,0.514,"cell recognition"),
c("GO:0046112","nucleobase biosynthetic process",0.099,0.3144,0.773,0.356,"cell recognition"),
c("GO:0006012","galactose metabolic process",0.121,0.4585,0.826,0.197,"cell recognition"),
c("GO:0006714","sesquiterpenoid metabolic process",0.181,0.9743,0.801,0.692,"cell recognition"),
c("GO:0006720","isoprenoid metabolic process",0.850,0.3105,0.786,0.549,"cell recognition"),
c("GO:0042214","terpene metabolic process",0.065,0.9831,0.813,0.168,"cell recognition"),
c("GO:0009069","serine family amino acid metabolic process",0.276,0.3357,0.771,0.668,"cell recognition"),
c("GO:0009072","aromatic amino acid family metabolic process",0.302,0.3597,0.755,0.387,"cell recognition"),
c("GO:0071554","cell wall organization or biogenesis",3.168,4.2151,0.919,0.061,"cell wall organization or biogenesis"),
c("GO:0005975","carbohydrate metabolic process",4.670,7.1916,0.877,0.062,"carbohydrate metabolism"),
c("GO:1901334","lactone metabolic process",0.026,1.5199,0.875,0.063,"lactone metabolism"),
c("GO:0009225","nucleotide-sugar metabolic process",0.207,0.9743,0.777,0.495,"lactone metabolism"),
c("GO:0006040","amino sugar metabolic process",0.138,0.4775,0.843,0.513,"lactone metabolism"),
c("GO:0016143","S-glycoside metabolic process",0.514,0.3019,0.746,0.597,"lactone metabolism"),
c("GO:1901136","carbohydrate derivative catabolic process",0.337,0.3707,0.808,0.552,"lactone metabolism"),
c("GO:0006022","aminoglycan metabolic process",0.086,0.9743,0.798,0.314,"lactone metabolism"),
c("GO:0042726","flavin-containing compound metabolic process",0.069,1.2312,0.837,0.146,"lactone metabolism"),
c("GO:0052646","alditol phosphate metabolic process",0.026,0.8575,0.831,0.428,"lactone metabolism"),
c("GO:0072593","reactive oxygen species metabolic process",0.734,0.4998,0.881,0.064,"reactive oxygen species metabolism"),
c("GO:0006793","phosphorus metabolic process",10.142,5.7295,0.849,0.114,"reactive oxygen species metabolism"),
c("GO:0018130","heterocycle biosynthetic process",14.014,0.7232,0.792,0.582,"reactive oxygen species metabolism"),
c("GO:0019438","aromatic compound biosynthetic process",14.247,0.9255,0.795,0.274,"reactive oxygen species metabolism"),
c("GO:1901362","organic cyclic compound biosynthetic process",14.769,0.4667,0.809,0.552,"reactive oxygen species metabolism"),
c("GO:0016310","phosphorylation",6.966,2.2582,0.846,0.166,"reactive oxygen species metabolism"),
c("GO:0016070","RNA metabolic process",17.562,0.4779,0.736,0.612,"reactive oxygen species metabolism"),
c("GO:0044260","cellular macromolecule metabolic process",35.028,3.6665,0.799,0.246,"reactive oxygen species metabolism"),
c("GO:0043412","macromolecule modification",15.900,0.3105,0.849,0.415,"reactive oxygen species metabolism"),
c("GO:0007154","cell communication",9.698,0.4576,0.913,0.076,"cell communication"),
c("GO:0001763","morphogenesis of a branching structure",0.125,0.9261,0.914,0.080,"morphogenesis of a branching structure"),
c("GO:0090603","sieve element differentiation",0.013,0.5291,0.881,0.237,"morphogenesis of a branching structure"),
c("GO:0048729","tissue morphogenesis",0.013,0.3597,0.920,0.414,"morphogenesis of a branching structure"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="autumn-winter-BP-revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
