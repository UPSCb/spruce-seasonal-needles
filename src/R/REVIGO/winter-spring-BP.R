

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
revigo.data <- rbind(c("GO:0005975","carbohydrate metabolic process",4.670,10.7378,0.856,0.000,"carbohydrate metabolism"),
c("GO:0050794","regulation of cellular process",22.244,0.8598,0.834,0.612,"carbohydrate metabolism"),
c("GO:0019538","protein metabolic process",19.012,0.9471,0.802,0.442,"carbohydrate metabolism"),
c("GO:0080090","regulation of primary metabolic process",13.681,1.1229,0.779,0.247,"carbohydrate metabolism"),
c("GO:0055065","metal ion homeostasis",0.842,0.7365,0.914,0.527,"carbohydrate metabolism"),
c("GO:0006793","phosphorus metabolic process",10.142,12.9841,0.837,0.140,"carbohydrate metabolism"),
c("GO:0006270","DNA replication initiation",0.086,0.9560,0.784,0.280,"carbohydrate metabolism"),
c("GO:0042967","acyl-carrier-protein biosynthetic process",2.838,1.3796,0.793,0.108,"carbohydrate metabolism"),
c("GO:0016070","RNA metabolic process",17.562,0.7617,0.706,0.430,"carbohydrate metabolism"),
c("GO:0044260","cellular macromolecule metabolic process",35.028,4.6224,0.766,0.260,"carbohydrate metabolism"),
c("GO:0043412","macromolecule modification",15.900,1.3647,0.813,0.415,"carbohydrate metabolism"),
c("GO:0010143","cutin biosynthetic process",0.091,0.5402,0.856,0.278,"carbohydrate metabolism"),
c("GO:0042446","hormone biosynthetic process",0.514,0.7735,0.837,0.241,"carbohydrate metabolism"),
c("GO:0010035","response to inorganic substance",3.699,3.7386,0.895,0.000,"response to inorganic substance"),
c("GO:0009628","response to abiotic stimulus",7.946,0.8058,0.907,0.433,"response to inorganic substance"),
c("GO:0001101","response to acid chemical",5.011,2.4088,0.893,0.638,"response to inorganic substance"),
c("GO:1901700","response to oxygen-containing compound",6.504,2.4088,0.891,0.605,"response to inorganic substance"),
c("GO:0009607","response to biotic stimulus",5.158,0.4583,0.910,0.396,"response to inorganic substance"),
c("GO:0042221","response to chemical",12.434,2.1046,0.905,0.371,"response to inorganic substance"),
c("GO:0015849","organic acid transport",0.682,5.1477,0.830,0.000,"organic acid transport"),
c("GO:0055085","transmembrane transport",5.175,2.8895,0.862,0.619,"organic acid transport"),
c("GO:0006811","ion transport",4.174,1.0691,0.863,0.519,"organic acid transport"),
c("GO:0030001","metal ion transport",1.588,1.5386,0.863,0.416,"organic acid transport"),
c("GO:0015696","ammonium transport",0.052,1.1229,0.876,0.611,"organic acid transport"),
c("GO:0033365","protein localization to organelle",1.148,0.9567,0.865,0.426,"organic acid transport"),
c("GO:0006857","oligopeptide transport",0.155,1.0459,0.874,0.476,"organic acid transport"),
c("GO:0071554","cell wall organization or biogenesis",3.168,4.3673,0.918,0.000,"cell wall organization or biogenesis"),
c("GO:0042026","protein refolding",0.039,0.4869,0.934,0.039,"protein refolding"),
c("GO:0007010","cytoskeleton organization",1.032,3.5018,0.892,0.054,"cytoskeleton organization"),
c("GO:0097435","supramolecular fiber organization",0.587,0.9471,0.830,0.361,"cytoskeleton organization"),
c("GO:0031163","metallo-sulfur cluster assembly",0.142,0.8598,0.903,0.315,"cytoskeleton organization"),
c("GO:0045229","external encapsulating structure organization",2.577,3.2013,0.886,0.426,"cytoskeleton organization"),
c("GO:0097264","self proteolysis",0.004,0.4869,0.896,0.062,"self proteolysis"),
c("GO:0006790","sulfur compound metabolic process",1.420,0.3547,0.864,0.069,"sulfur compound metabolism"),
c("GO:0001763","morphogenesis of a branching structure",0.125,0.6409,0.926,0.087,"morphogenesis of a branching structure"),
c("GO:0009225","nucleotide-sugar metabolic process",0.207,4.0068,0.797,0.088,"nucleotide-sugar metabolism"),
c("GO:0033865","nucleoside bisphosphate metabolic process",0.069,0.7078,0.741,0.430,"nucleotide-sugar metabolism"),
c("GO:0055086","nucleobase-containing small molecule metabolic process",1.800,1.5991,0.693,0.236,"nucleotide-sugar metabolism"),
c("GO:0042816","vitamin B6 metabolic process",0.039,0.3438,0.761,0.409,"nucleotide-sugar metabolism"),
c("GO:0046112","nucleobase biosynthetic process",0.099,0.3334,0.721,0.556,"nucleotide-sugar metabolism"),
c("GO:1901334","lactone metabolic process",0.026,2.1046,0.858,0.157,"nucleotide-sugar metabolism"),
c("GO:0006022","aminoglycan metabolic process",0.086,1.4801,0.795,0.495,"nucleotide-sugar metabolism"),
c("GO:0030243","cellulose metabolic process",0.337,1.4021,0.807,0.471,"nucleotide-sugar metabolism"),
c("GO:0019694","alkanesulfonate metabolic process",4.265,0.9571,0.771,0.385,"nucleotide-sugar metabolism"),
c("GO:0006012","galactose metabolic process",0.121,1.5991,0.753,0.452,"nucleotide-sugar metabolism"),
c("GO:0046184","aldehyde biosynthetic process",0.052,0.7533,0.799,0.239,"nucleotide-sugar metabolism"),
c("GO:0006714","sesquiterpenoid metabolic process",0.181,1.4598,0.788,0.268,"nucleotide-sugar metabolism"),
c("GO:0042726","flavin-containing compound metabolic process",0.069,0.4841,0.811,0.168,"nucleotide-sugar metabolism"),
c("GO:0052646","alditol phosphate metabolic process",0.026,0.3819,0.853,0.454,"nucleotide-sugar metabolism"),
c("GO:0055114","oxidation-reduction process",7.462,1.1229,0.777,0.481,"nucleotide-sugar metabolism"),
c("GO:0042214","terpene metabolic process",0.065,0.7629,0.801,0.692,"nucleotide-sugar metabolism"),
c("GO:0009069","serine family amino acid metabolic process",0.276,0.3467,0.725,0.668,"nucleotide-sugar metabolism"),
c("GO:0009072","aromatic amino acid family metabolic process",0.302,0.4869,0.717,0.494,"nucleotide-sugar metabolism"),
c("GO:0019499","cyanide metabolic process",0.009,0.9471,0.819,0.132,"nucleotide-sugar metabolism"),
c("GO:0044262","cellular carbohydrate metabolic process",1.817,0.4869,0.800,0.602,"nucleotide-sugar metabolism"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="winter-spring-BP-revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

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
