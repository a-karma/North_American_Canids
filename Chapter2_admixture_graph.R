#-------------------------------------------------------------------LOADING-----------------------------------------------------------------
# install.packages("admixturegraph", dependencies = T)
setwd("Documents/Canids/2017-10-15 Admixture project/Admixture Graph") 
getwd()

# loading required package
library(admixturegraph)

############################################################## Loading data (outgroup: andean fox)
d_stats_ch4 <- read.table(file = "d_stats/ch4_dstat.txt", header =T)
d_stats_all_ch <- read.table(file = "d_stats/all_ch_dstat.txt", header =T)
d_stats_autosomes <- read.table(file = "d_stats/autosomes_dstat.txt", header =T)
f4_autosomes<- read.table(file = "d_stats/autosomes_f4.txt", header =T)
# defining LEAVES
outgroup <- "Fox"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

############################################################### Loading data (outgroup: jackal, American wolf: yellowstone)
f4_ch1_jack <- read.table(file = "d_stats/jack_ch1_f4.txt", header =T)
f4_auto_jack <- read.table(file = "d_stats/autosomes_f4_jack.txt", header =T)
# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

################################################################ Loading data (outgroup: jackal, American wolf: Banks Island)
f4_Banks_auto <- read.table(file = "d_stats/auto_Banks_f4.txt", header =T)
f4_30kb_Banks <- read.table(file = "d_stats/30kb_Banks_f4.txt", header =T)
# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

################################################################# Loading data (outgroup jackal and pruned snps dataset)
f4_100kb <- read.table(file = "d_stats/100kb_spaced_f4.txt", header =T)
f4_50kb <- read.table(file = "d_stats/50kb_spaced_f4.txt", header =T)
f4_30kb <- read.table(file = "d_stats/30kb_spaced_f4.txt", header =T)
f4_10kb <- read.table(file = "d_stats/10kb_spaced_f4.txt", header =T)
# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

################################################################# Different reference genome (African Wild Dog)
f4_AWD_5 <- read.table(file = "d_stats/AWD_5_f4.txt", header =T)
f4_AWD_30kb <- read.table(file = "d_stats/30kb_awd_f4.txt", header =T)
# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

################################################################ Different ref gen (Red Fox)
f4_RF_5 <- read.table(file = "d_stats/RF_5_f4.txt", header =T)
f4_RF_30kb <- read.table(file = "d_stats/RF_5_30kb_f4.txt", header =T)
f4_RF_5_filt <- read.table(file = "d_stats/RF_5_filtered_f4.txt", header =T)
f4_RF_filt_30kb <- read.table(file = "d_stats/RF_5_filtered_30kb_f4.txt", header =T)
# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

################################################################ Loading data (outgroup jack and increasing filtering threshold)
f4_fil2 <- read.table(file = "d_stats/filter_ex_cpg_f4.txt", header =T)
f4_fil3 <- read.table(file = "d_stats/filter_ex_cpg_rep_f4.txt", header =T)
f4_fil4 <- read.table(file = "d_stats/filter_ex_cpg_rep_low_f4.txt", header =T)
f4_fil2_auto <- read.table(file = "d_stats/filter_ex_cpg_auto_f4.txt", header =T)
f4_fil3_auto <- read.table(file = "d_stats/filter_ex_cpg_rep_auto_f4.txt", header =T)
f4_fil4_auto <- read.table(file = "d_stats/filter_ex_cpg_rep_low_auto_f4.txt", header =T)
f4_fil2_auto_30kb_sp <- read.table(file = "d_stats/filter_ex_cpg_auto_30kb_f4.txt", header =T)
f4_fil3_auto_30kb_sp <- read.table(file = "d_stats/filter_ex_cpg_rep_auto_30kb_f4.txt", header =T)

# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote","Red", outgroup)

#################################################################### Variable Assignment
canids<-f4_30kb_Banks

#---------------------------------------------------------------Divergence Graph--------------------------------------------------------------
# variables names are pretty self-explanatory. Note the sintax for the edge construction: backward in time.
# i.e. the first element must be the node closer to the present (child) and the second element must be the inner node (parent).

# Divergence model
inner_nodes <- c("R","Holo_W","WC","RC")
edges <- parent_edges(c(edge(outgroup,"R"),
                        edge("WC","R"),
                        edge("E_wolf","Holo_W"),
                        edge("A_wolf","Holo_W"),
                        edge("Holo_W","WC"),
                        edge("Coyote","RC"),
                        edge("Red","RC"),
                        edge("RC","WC")))
null_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the null graph to our data
null_k9_fit <- fit_graph(canids, null_k9)
summary(null_k9_fit)
# plotting results
plot(null_k9, show_inner_node_labels = TRUE)
plot(null_k9_fit)

# Alternative Divergence Model
inner_nodes <- c("R","Holo_W","WC","WR")
edges_alt <- parent_edges(c(edge(outgroup,"R"),
                        edge("WC","R"),
                        edge("E_wolf","Holo_W"),
                        edge("A_wolf","Holo_W"),
                        edge("Holo_W","WR"),
                        edge("Coyote","WC"),
                        edge("Red","WR"),
                        edge("WR","WC")))

null_alt_k9 <- agraph(leaves, inner_nodes, edges_alt)
# fitting the alternative null graph to our data
null_alt_k9_fit <- fit_graph(canids, null_alt_k9)
summary(null_alt_k9_fit)
# plotting results
plot(null_alt_k9, show_inner_node_labels = TRUE)
plot(null_alt_k9_fit)

# Polytomy model
inner_nodes <- c("R","Holo_W","WC")
edges_poly <- parent_edges(c(edge(outgroup,"R"),
                            edge("WC","R"),
                            edge("E_wolf","Holo_W"),
                            edge("A_wolf","Holo_W"),
                            edge("Holo_W","WC"),
                            edge("Coyote","WC"),
                            edge("Red","Holo_W")))
null_poly_k9 <- agraph(leaves, inner_nodes, edges_poly)
# fitting the alternative null graph to our data
null_poly_k9_fit <- fit_graph(canids, null_poly_k9)
summary(null_poly_k9_fit)
# plotting results
plot(null_poly_k9, show_inner_node_labels = TRUE)
plot(null_poly_k9_fit)

# Very recent divergence model
inner_nodes <- c("R","Holo_W","WC","AAW")
edges_rec <- parent_edges(c(edge(outgroup,"R"),
                            edge("WC","R"),
                            edge("E_wolf","Holo_W"),
                            edge("AAW","Holo_W"),
                            edge("A_wolf","AAW"),
                            edge("Holo_W","WC"),
                            edge("Coyote","WC"),
                            edge("Red","AAW")))
null_recent_k9 <- agraph(leaves, inner_nodes, edges_rec)
# fitting the alternative null graph to our data
null_recent_k9_fit <- fit_graph(canids, null_recent_k9)
summary(null_recent_k9_fit)
# plotting results
plot(null_recent_k9, show_inner_node_labels = TRUE)
plot(null_recent_k9_fit)


#------------------------------------------------------------Single Admixture Event-----------------------------------------------------------
# including admixtures: same leaves but different edges and inner nodes

# Recent admiture model (two null models nested in here: null_recent_k9 and null_k9, intuitively mergin the two)
inner_nodes <- c("R","Holo_W","WC","AAW","RA","RC")
edges_rec_adm <- parent_edges(c(edge(outgroup,"R"),
                            edge("WC","R"),
                            edge("E_wolf","Holo_W"),
                            edge("AAW","Holo_W"),
                            edge("A_wolf","AAW"),
                            edge("Holo_W","WC"),
                            edge("Coyote","RC"),
                            edge("RC","WC"),
                            edge("Red","RA"),
                            admixture_edge("RA","RC","AAW","alpha")))
rec_adm_k9 <- agraph(leaves, inner_nodes, edges_rec_adm)
# fitting the alternative null graph to our data
rec_adm_k9_fit <- fit_graph(canids, rec_adm_k9)
summary(rec_adm_k9_fit)
# plotting results
plot(rec_adm_k9, show_inner_node_labels = TRUE)
plot(rec_adm_k9_fit)

# older admixture event (which is exactly the same as adding an admixture edge to the alternative div graph)
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW")
edges_adm <- parent_edges(c(edge("E_wolf","Holo_W"),
                            edge("A_wolf","Holo_W"),
                            edge("Holo_W","AHW"),
                            edge("AHW","WC"),
                            edge("Coyote","RC"),
                            edge("Red","RA"),
                            edge("RC","WC"),
                            edge("WC","Root"),
                            edge(outgroup,"Root"),
                            admixture_edge("RA","RC","AHW","alpha")))
old_adm_k9 <- agraph(leaves, inner_nodes, edges_adm)
# model fitting
old_adm_k9_fit <- fit_graph(canids, old_adm_k9)
summary(old_adm_k9_fit)
# plotting results
plot(old_adm_k9, show_inner_node_labels = TRUE)
plot(old_adm_k9_fit)

# variant older admixture event (which is exactly the same as adding an admixture edge to the recent div graph)
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","ghost")
edges_adm_r <- parent_edges(c(edge("E_wolf","Holo_W"),
                            edge("A_wolf","ghost"),
                            edge("Holo_W","AHW"),
                            edge("AHW","WC"),
                            edge("Coyote","RC"),
                            edge("ghost","RA"),
                            edge("Red","ghost"),
                            edge("RC","WC"),
                            edge("WC","Root"),
                            edge(outgroup,"Root"),
                            admixture_edge("RA","RC","AHW","alpha")))
old_adm_var_k9 <- agraph(leaves, inner_nodes, edges_adm_r)
# model fitting
old_adm_var_k9_fit <- fit_graph(canids, old_adm_var_k9)
summary(old_adm_var_k9_fit)
# plotting results
plot(old_adm_var_k9, show_inner_node_labels = TRUE)
plot(old_adm_var_k9_fit)

#----------------------------------------------------------------Double Admixture-------------------------------------------------------------
# old admixture plus recen backflow into american wolf
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA")
edges_adm2 <- parent_edges(c(edge("E_wolf","Holo_W"),
                            edge("A_wolf","AAW"),
                            edge("Holo_W","AHW"),
                            edge("AHW","WC"),
                            edge("Coyote","RC"),
                            edge("RA","RAA"),
                            edge("Red","RA"),
                            edge("RC","WC"),
                            edge("WC","Root"),
                            edge(outgroup,"Root"),
                            admixture_edge("RAA","RC","AHW","alpha"),
                            admixture_edge("AAW","RA","Holo_W","beta")))
double_adm_k9 <- agraph(leaves, inner_nodes, edges_adm2)
# model fitting
double_adm_k9_fit <- fit_graph(canids, double_adm_k9)
summary(double_adm_k9_fit)
# plotting results
plot(double_adm_k9)
plot(double_adm_k9_fit)

# coyote introgression into American Wolf
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","CA")
edges_adm_cia <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","WC"),
                             edge("Coyote","CA"),
                             edge("CA","RC"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RA","RC","AHW","alpha"),
                             admixture_edge("AAW","CA","Holo_W","beta")))
double_adm_cia_k9 <- agraph(leaves, inner_nodes, edges_adm_cia)
# model fitting
double_adm_cia_k9_fit <- fit_graph(canids, double_adm_cia_k9)
summary(double_adm_cia_k9_fit)
# plotting results
plot(double_adm_cia_k9,show_inner_node_labels = TRUE)
plot(double_adm_cia_k9_fit)

# old admixture and recent coyote introgression into the red wolf lineage
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","RAA","CA")
edges_adm_cir <- parent_edges(c(edge("E_wolf","Holo_W"),
                                edge("A_wolf","Holo_W"),
                                edge("Holo_W","AHW"),
                                edge("AHW","WC"),
                                edge("Coyote","CA"),
                                edge("CA","RC"),
                                edge("Red","RA"),
                                edge("RC","WC"),
                                edge("WC","Root"),
                                edge(outgroup,"Root"),
                                admixture_edge("RAA","RC","AHW","alpha"),
                                admixture_edge("RA","RAA","CA","beta")))
double_adm_cir_k9 <- agraph(leaves, inner_nodes, edges_adm_cir)
# model fitting
double_adm_cir_k9_fit <- fit_graph(canids, double_adm_cir_k9)
summary(double_adm_cir_k9_fit)
# plotting results
plot(double_adm_cir_k9,show_inner_node_labels = TRUE)
plot(double_adm_cir_k9_fit)

# old admixture and recent introgression of American wolves into the red wolf lineage
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA")
edges_adm2b <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("A_wolf","AAW"),
                             edge("AAW","Holo_W"),
                             edge("Holo_W","AHW"),
                             edge("AHW","WC"),
                             edge("Coyote","RC"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("RA","RAA","AAW","beta")))
double_adm2_k9 <- agraph(leaves, inner_nodes, edges_adm2b)
# model fitting
double_adm2_k9_fit <- fit_graph(canids, double_adm2_k9)
summary(double_adm2_k9_fit)
# plotting results
plot(double_adm2_k9,show_inner_node_labels = TRUE)
plot(double_adm2_k9_fit)

#par(mfcol=c(3,2)) # to plot non-fitting models
#plot(null_k9, show_inner_node_labels = TRUE)
#plot(null_alt_k9, show_inner_node_labels = TRUE)
#plot(null_poly_k9, show_inner_node_labels = TRUE)
#plot(null_recent_k9, show_inner_node_labels = TRUE)
#plot(rec_adm_k9, show_inner_node_labels = TRUE)
#plot(old_adm_k9, show_inner_node_labels = TRUE)
#plot(old_adm_var_k9, show_inner_node_labels = TRUE)
#plot(double_adm_cia_k9,show_inner_node_labels = TRUE)
#plot(double_adm_cir_k9,show_inner_node_labels = TRUE)
#plot(double_adm2_k9,show_inner_node_labels = TRUE)


#------------------------------------------------------------Three Admixture Events----------------------------------------------------------
# old admixture plus recen backflow into american wolf and of coyote into the red
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA","CA","Hist_R")
edges_adm3 <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","WC"),
                             edge("CA","RC"),
                             edge("Coyote","CA"),
                             edge("RA","RAA"),
                             edge("Red","Hist_R"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta"),
                             admixture_edge("Hist_R","RA","CA")))
tri_adm_k9 <- agraph(leaves, inner_nodes, edges_adm3)
# model fitting
tri_adm_k9_fit <- fit_graph(canids, tri_adm_k9)
summary(tri_adm_k9_fit)
# plotting results
plot(tri_adm_k9,show_inner_node_labels = TRUE)
plot(tri_adm_k9_fit)

# old admixture plus recen backflow into american wolf followe by wolves geen flow
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA","Ghost1","Hist_AW")
edges_adm_3crazy <- parent_edges(c(edge("E_wolf","Ghost1"),
                             edge("Ghost1","Holo_W"),
                             edge("A_wolf","Hist_AW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","WC"),
                             edge("Coyote","RC"),
                             edge("RA","RAA"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta"),
                             admixture_edge("Hist_AW","Ghost1","AAW","gamma")))
crazy3_adm_k9 <- agraph(leaves, inner_nodes, edges_adm_3crazy)
crazy3_adm_k9_fit <- fit_graph(canids, crazy3_adm_k9)
summary(crazy3_adm_k9_fit)
# plotting results
plot(crazy3_adm_k9,show_inner_node_labels = TRUE)
plot(crazy3_adm_k9_fit)


###############################################################################################################################################

#-----------------------------------------------------------Including Ancient Genome-----------------------------------------------------------

# Full Dataset
f4_6_samples <- read.table(file = "d_stats/6_samples_f4.txt", header =T)
f4_6_samples_v2 <-read.table(file = "d_stats/6_samples_v2_f4.txt", header =T)

# Pruned Dataset
f4_6_pruned_10kb <- read.table(file = "d_stats/6_pruned_10kb_f4.txt", header =T)
f4_6_pruned_n50 <- read.table(file = "d_stats/6_pruned_f4.txt", header =T)
f4_6_30kb_spaced <- read.table(file = "d_stats/6_30kb_spaced_f4.txt", header =T)

# defining LEAVES
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Ancient","Coyote","Red", outgroup)


canids <- f4_6_30kb_spaced

# Divergence model
inner_nodes <- c("R","Holo_W","Pleisto_W","WC","RC")
edges <- parent_edges(c(edge(outgroup,"R"),
                        edge("Ancient","Pleisto_W"),
                        edge("WC","R"),
                        edge("E_wolf","Holo_W"),
                        edge("A_wolf","Holo_W"),
                        edge("Holo_W","Pleisto_W"),
                        edge("Pleisto_W","WC"),
                        edge("Coyote","RC"),
                        edge("Red","RC"),
                        edge("RC","WC")))
null_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the null graph to our data
null_k9_fit <- fit_graph(canids, null_k9)
summary(null_k9_fit)
# plotting results
plot(null_k9, show_inner_node_labels = TRUE)
plot(null_k9_fit)


# Alternative Divergence model
inner_nodes <- c("R","Holo_W","Pleisto_W","WC","RC")
edges <- parent_edges(c(edge(outgroup,"R"),
                        edge("Ancient","Holo_W"),
                        edge("WC","R"),
                        edge("E_wolf","Holo_W"),
                        edge("A_wolf","Pleisto_W"),
                        edge("Holo_W","Pleisto_W"),
                        edge("Pleisto_W","WC"),
                        edge("Coyote","RC"),
                        edge("Red","RC"),
                        edge("RC","WC")))
null_alt_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the null graph to our data
null_alt_k9_fit <- fit_graph(canids, null_alt_k9)
summary(null_k9_fit)
# plotting results
plot(null_alt_k9, show_inner_node_labels = TRUE)
plot(null_alt_k9_fit)


# old admixture plus recen backflow into american wolf (three version of this model)
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA","Pleisto_W")
edges_adm2 <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("Ancient","Pleisto_W"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","Pleisto_W"),
                             edge("Pleisto_W","WC"),
                             edge("Coyote","RC"),
                             edge("RA","RAA"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta")))
double_adm_k9 <- agraph(leaves, inner_nodes, edges_adm2)
# model fitting
double_adm_k9_fit <- fit_graph(canids, double_adm_k9)
summary(double_adm_k9_fit)
# plotting results
plot(double_adm_k9,show_inner_node_labels = TRUE)
plot(double_adm_k9_fit)

inner_nodes <- c("Root","Eur_Pleisto","WC","RC","RA","Anc_Am_Pleisto","AAW","RAA","A_Pleisto_W")
edges_adm2_v1 <- parent_edges(c(edge("E_wolf","Eur_Pleisto"),
                                edge("Ancient","Eur_Pleisto"),
                                edge("A_wolf","AAW"),
                                edge("Eur_Pleisto","Anc_Am_Pleisto"),
                                edge("Anc_Am_Pleisto","A_Pleisto_W"),
                                edge("A_Pleisto_W","WC"),
                                edge("Coyote","RC"),
                                edge("RA","RAA"),
                                edge("Red","RA"),
                                edge("RC","WC"),
                                edge("WC","Root"),
                                edge(outgroup,"Root"),
                                admixture_edge("RAA","RC","A_Pleisto_W","alpha"),
                                admixture_edge("AAW","RA","Anc_Am_Pleisto","beta")))
double_adm_v1_k9 <- agraph(leaves, inner_nodes, edges_adm2_v1)
# model fitting
double_adm_v1_k9_fit <- fit_graph(canids, double_adm_v1_k9)
summary(double_adm_v1_k9_fit)
# plotting results
plot(double_adm_v1_k9,show_inner_node_labels = TRUE)
plot(double_adm_v1_k9_fit)

inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA","Pleisto_W")
edges_adm2_v2 <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("Ancient","AHW"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","Pleisto_W"),
                             edge("Pleisto_W","WC"),
                             edge("Coyote","RC"),
                             edge("RA","RAA"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","Pleisto_W","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta")))
double_adm_v2_k9 <- agraph(leaves, inner_nodes, edges_adm2_v2)
# model fitting
double_adm_v2_k9_fit <- fit_graph(canids, double_adm_v2_k9)
summary(double_adm_v2_k9_fit)
# plotting results
plot(double_adm_v2_k9,show_inner_node_labels = TRUE)
plot(double_adm_v2_k9_fit)


# three admixture: pleistocene non-sense
inner_nodes <- c("Root","Eur_Pleisto","holo","Am_pleisto","WC","RC","RA","Pleistocene wolves","AAW","RAA","MRCA wolves")
edges_adm3_v2 <- parent_edges(c(edge("E_wolf","Eur_Pleisto"),
                                edge("A_wolf","AAW"),
                                edge("Ancient","holo"),
                                edge("Am_pleisto","Pleistocene wolves"),
                                edge("Eur_Pleisto","Pleistocene wolves"),
                                edge("Pleistocene wolves","MRCA wolves"),
                                edge("MRCA wolves","WC"),
                                edge("Coyote","RC"),
                                edge("RA","RAA"),
                                edge("Red","RA"),
                                edge("RC","WC"),
                                edge("WC","Root"),
                                edge(outgroup,"Root"),
                                admixture_edge("RAA","RC","MRCA wolves","alpha"),
                                admixture_edge("AAW","RA","Am_pleisto","beta"),
                                admixture_edge("holo","Am_pleisto","Eur_Pleisto","gamma")))
triple_adm_v2_k9 <- agraph(leaves, inner_nodes, edges_adm3_v2)
# model fitting
triple_adm_v2_k9_fit <- fit_graph(canids, triple_adm_v2_k9)
summary(triple_adm_v2_k9_fit)
# plotting results
plot(triple_adm_v2_k9,show_inner_node_labels = TRUE)
plot(triple_adm_v2_k9_fit)

# three admixture: old admixture, pleistocene loop, and recen backflow into american wolf
inner_nodes <- c("Root","Holo_W","WC","RC","RA","AHW","AAW","RAA","MRCA wolves","Pleisto_W","late_pleisto")
edges_adm3l <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("Ancient","Pleisto_W"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","late_pleisto"),
                             edge("Pleisto_W","MRCA wolves"),
                             edge("AHW","MRCA wolves"),
                             edge("MRCA wolves","WC"),
                             edge("Coyote","RC"),
                             edge("RA","RAA"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta"),
                             admixture_edge("late_pleisto","AHW","Pleisto_W","gamma")))
triple_adm_loop_k9 <- agraph(leaves, inner_nodes, edges_adm3l)
# model fitting
triple_adm_loop_k9_fit <- fit_graph(canids, triple_adm_loop_k9)
summary(triple_adm_loop_k9_fit)
# plotting results
plot(triple_adm_loop_k9, show_inner_node_labels = TRUE)
plot(triple_adm_loop_k9_fit)


#par(mfcol=c(3,2))
#plot(null_k9, show_inner_node_labels = TRUE)
#plot(null_alt_k9, show_inner_node_labels = TRUE)
#plot(double_adm_k9,show_inner_node_labels = TRUE)
#plot(double_adm_v1_k9,show_inner_node_labels = TRUE)
#plot(triple_adm_v2_k9,show_inner_node_labels = TRUE)
#plot(triple_adm_loop_k9, show_inner_node_labels = TRUE)
#par()

#############################################################################################################################################

#--------------------------------------------------Including Ancient and Historical Genomes--------------------------------------------------

f4_7_samples <- read.table(file = "d_stats/7_jack_f4.txt", header =T)
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Ancient","Red","Coyote","ancoy1", outgroup)
canids <- f4_7_samples

# Divergence Model
inner_nodes <- c("R","Holo_W","Pleisto_W","AC","WC","RC")
edges <- parent_edges(c(edge(outgroup,"R"),
                        edge("Ancient","Pleisto_W"),
                        edge("WC","R"),
                        edge("E_wolf","Holo_W"),
                        edge("A_wolf","Holo_W"),
                        edge("Holo_W","Pleisto_W"),
                        edge("Pleisto_W","WC"),
                        edge("Coyote","AC"),
                        edge("ancoy1","AC"),
                        edge("AC","RC"),
                        edge("Red","RC"),
                        edge("RC","WC")))
null_k9_7 <- agraph(leaves, inner_nodes, edges)
#fitting the null graph to our data
null_k9_7_fit <- fit_graph(canids, null_k9_7)
summary(null_k9_7_fit)
# plotting results
plot(null_k9_7, show_inner_node_labels = TRUE)
plot(null_k9_7_fit)

# old admixture plus recen backflow into american wolf (three version of this model)
inner_nodes <- c("Root","Holo_W","WC","AC","RC","RA","AHW","AAW","RAA","Pleisto_W")
edges_adm2_7 <- parent_edges(c(edge("E_wolf","Holo_W"),
                             edge("Ancient","Pleisto_W"),
                             edge("A_wolf","AAW"),
                             edge("Holo_W","AHW"),
                             edge("AHW","Pleisto_W"),
                             edge("Pleisto_W","WC"),
                             edge("Coyote","AC"),
                             edge("ancoy1","AC"),
                             edge("AC","RC"),
                             edge("RA","RAA"),
                             edge("Red","RA"),
                             edge("RC","WC"),
                             edge("WC","Root"),
                             edge(outgroup,"Root"),
                             admixture_edge("RAA","RC","AHW","alpha"),
                             admixture_edge("AAW","RA","Holo_W","beta")))
double_adm_7_k9 <- agraph(leaves, inner_nodes, edges_adm2_7)
# model fitting
double_adm_7_k9_fit <- fit_graph(canids, double_adm_7_k9)
summary(double_adm_7_k9_fit)
# plotting results
plot(double_adm_7_k9,show_inner_node_labels = TRUE)
plot(double_adm_7_k9_fit)

#############################################################################################################################################

#-------------------------------------------------------------Poster Figures-----------------------------------------------------------------
outgroup <- "Jackal"
leaves <- c("E_wolf", "A_wolf","Coyote", outgroup)
inner_nodes <- c("Root","Wolves Ancestor","Wolf/Coyote Ancestor")
edges <- parent_edges(c(edge(outgroup,"Root"),
                        edge("Wolf/Coyote Ancestor","Root"),
                        edge("E_wolf","Wolves Ancestor"),
                        edge("A_wolf","Wolves Ancestor"),
                        edge("Wolves Ancestor","Wolf/Coyote Ancestor"),
                        edge("Coyote","Wolf/Coyote Ancestor")))
div_EAwolf <- agraph(leaves, inner_nodes, edges)
plot(div_EAwolf, show_inner_node_labels = TRUE,platform=0.5)


outgroup <- "Jackal"
leaves <- c("A_wolf","Red Wolf","Coyote", outgroup)
inner_nodes <- c("Root","Red/Coyote Ancestor","Wolf/Coyote Ancestor")
edges <- parent_edges(c(edge(outgroup,"Root"),
                        edge("Wolf/Coyote Ancestor","Root"),
                        edge("A_wolf","Wolf/Coyote Ancestor"),
                        edge("Red/Coyote Ancestor","Wolf/Coyote Ancestor"),
                        edge("Coyote","Red/Coyote Ancestor"),
                        edge("Red Wolf","Red/Coyote Ancestor")))
div_Red <- agraph(leaves, inner_nodes, edges)
plot(div_Red, show_inner_node_labels = TRUE)

##############################################################################################################################################

              #########
              # Notes #
              #########

# 1) Five sample dataset
# The best fitting model is always the double admixture model with ancient admixture plus recent backflow
# This pattern is consistent in both outgroup cases.
# The reference genome seems not to have an impact on this (tested both AWD and RedFox)
# Excluding the X chromosome doesn't have any effect (except reducing the total number of loci).
# The pattern looks consistent regardless the American wolf used.
# Providing that the pruning is done by spacing snps, it does not have a significant impact.
# Filtering out exon and cpg island slightly reduced the number of snps in the dataset but didn't alter the results.
# Adding another filtering step by removing all repetitive sequences also didn't alter the results.
# Even when combining these three filtering stages with pruning snps based on their physical distance and excluding the X, same best fit
# Filtering also region with low mappability removed most of the snps from the dataset reducing the power of the analysis and
# possibly artificially enrich the amount of snps in more conserved regions. 
# In this case, the Recent single Admiture Model, 
# the old admixture and recent coyote introgression into the red wolf lineage model, 
# the ancient admixture plus recent backflow from Red into American wolf model, 
# the old admixture and recent introgression of American wolves into the red wolf lineage
# are indistinguishable and they all fitt equally well.


# 2) Six sample dataset
# Including an ancient genome slightly complicated the analysis. The best fitting model is still the same although the fit is not perfect.
# This happens regardless the reference genome or the spacing















