library(dplyr)
library(tidyr)
library(lmerTest)
library(igraph)


# target_species <- c(
#   "s__Fusobacterium_polymorphum", 
#   "s__GGB49229_SGB69060", 
#   "s__Leptotrichia_massiliensis", 
#   "s__Leptotrichia_wadei", 
#   "s__Prevotella_veroralis", 
#   "s__Stomatobaculum_SGB5266", 
#   "s__Fusobacterium_pseudoperiodonticum", 
#   "s__Leptotrichia_sp_oral_taxon_212", 
#   "s__GGB49401_SGB69314", 
#   "s__GGB74353_SGB98242", 
#   "s__Lachnoanaerobaculum_sp_ICM7", 
#   "s__Leptotrichia_sp_oral_taxon_221", 
#   "s__Pseudoramibacter_SGB4082", 
#   "s__Prevotella_histicola"
# )

target_species <- c(
  "s__Peptidiphaga_SGB15895",
  "s__GGB6814_SGB9728",
  "s__Selenomonas_flueggei",
  "s__Lachnoanaerobaculum_orale",
  "s__GGB2666_SGB3592",
  "s__GGB3388_SGB4476",
  "s__GGB3059_SGB69431",
  "s__Catonella_morbi",
  "s__Lachnoanaerobaculum_sp_ICM7",
  "s__Neisseria_bacilliformis",
  "s__TM7_phylum_sp_oral_taxon_348",
  "s__Leptotrichia_massiliensis",
  "s__Leptotrichia_sp_oral_taxon_215",
  "s__Lachnoanaerobaculum_sp_Marseille_Q4761",
  "s__Leptotrichia_hongkongensis",
  "s__Fusobacterium_animalis",
  "s__Leptotrichia_hofstadii",
  "s__Leptotrichia_wadei",
  "s__Fusobacterium_polymorphum",
  "s__Capnocytophaga_genosp_AHN8471"
)



target_mets <- c("1-palmitoyl-sn-glycero-3-phosphocholine", "Lpc(18:1)", "Phosphatidylcholine lyso 18:0", "Pc(18:2/0:0)", "Lpc(20:1-sn1)", "Pc(17:0/0:0)")
target_clinicals <- c("PLI", "GI", "BI")

# =========================================================

abd_df <- as.data.frame(t(read.table("Microbial_Relative_Abundance_Matrix.txt", header=T, row.names=1, sep="\t", check.names=F)), check.names=F) %>% 
  mutate(SampleID = rownames(.)) %>% inner_join(read.csv("metadata2.csv", check.names=F) %>% select(SampleID, Subject_ID, Time), by="SampleID") %>% select(-SampleID)
colnames(abd_df) <- sub(".*\\|(s__.*)", "\\1", colnames(abd_df))


met_df <- as.data.frame(t(read.csv("blood_metabolitis2.csv", row.names=1, check.names=F)), check.names=F) %>% 
  mutate(Sample_ID = rownames(.), Subject_ID = sub("_.*", "", Sample_ID), TimeCode = sub(".*_", "", Sample_ID), Time = case_when(TimeCode=="1"~"W0", TimeCode=="2"~"W2", TimeCode=="3"~"W3")) %>% select(-Sample_ID, -TimeCode)

pli <- read.csv("data_PLI_2.csv", check.names=F) %>% pivot_longer(-Subject_ID, names_to="Time", values_to="PLI") %>% mutate(Time=sub(".*_", "", Time))
gi <- read.csv("data_GI_2.csv", check.names=F) %>% pivot_longer(-Subject_ID, names_to="Time", values_to="GI") %>% mutate(Time=sub(".*_", "", Time))
bi <- read.csv("data_BI_2.csv", check.names=F) %>% pivot_longer(-Subject_ID, names_to="Time", values_to="BI") %>% mutate(Time=sub(".*_", "", Time))
clin_df <- pli %>% full_join(gi, by=c("Subject_ID", "Time")) %>% full_join(bi, by=c("Subject_ID", "Time"))

meta.tab <- inner_join(read.csv("metadata__blood2.csv", check.names=F), inner_join(abd_df, met_df, by=c("Subject_ID", "Time")), by="Subject_ID") %>%
  inner_join(clin_df, by=c("Subject_ID", "Time"))
meta.tab$Group <- factor(meta.tab$Group, levels=c("Experimental", "Control"), labels=c("S", "C"))
meta.tab$Age <- as.numeric(meta.tab$Age)


for (var in c(target_species, target_clinicals)) {
  if(var %in% colnames(meta.tab)) {
    val <- as.numeric(meta.tab[[var]])
    min_v <- min(val[val>0], na.rm=T) / 2
    if(is.infinite(min_v) || is.na(min_v)) min_v <- 1e-6
    meta.tab[[paste0("z_", var)]] <- as.numeric(scale(log2(val + min_v)))
  }
}


for (var in target_mets) {
  if(var %in% colnames(meta.tab)) {
    meta.tab[[paste0("z_", var)]] <- as.numeric(scale(as.numeric(meta.tab[[var]])))
  }
}


# =========================================================

get_network_edges <- function(var_list1, var_list2, group_name) {
  edges <- data.frame()
  df_sub <- meta.tab %>% filter(Group == group_name)
  
  for(v1 in var_list1) {
    for(v2 in var_list2) {
      if(!paste0("z_",v1) %in% colnames(meta.tab) || !paste0("z_",v2) %in% colnames(meta.tab)) next
      
      f_grp <- tryCatch(lmerTest::lmer(as.formula(paste0("`z_", v1, "` ~ `z_", v2, "` + Age + Gender + (1|Subject_ID)")), data = df_sub), error=function(e) NULL)
      grp_p <- NA; grp_coef <- NA
      if(!is.null(f_grp)) { 
        cf <- summary(f_grp)$coefficients
        term <- grep(paste0("z_",v2), rownames(cf), value=T, fixed=T)[1]
        if(!is.na(term)) { grp_p <- cf[term, "Pr(>|t|)"]; grp_coef <- cf[term, "Estimate"] } 
      }
      
      edges <- rbind(edges, data.frame(node1=v1, node2=v2, coef=grp_coef, pval=grp_p))
    }
  }
  
  if(nrow(edges) > 0) {
    edges$qval <- p.adjust(edges$pval, method="BH") 
  }
  return(edges)
}

# =========================================================

export_network_for_group <- function(group_name) {

  
  edge_spe_met <- get_network_edges(target_species, target_mets, group_name)
  edge_met_clin <- get_network_edges(target_mets, target_clinicals, group_name)
  edge_spe_clin <- get_network_edges(target_species, target_clinicals, group_name)
  

  all_edges <- rbind(edge_spe_met, edge_met_clin, edge_spe_clin) %>%
    filter(!is.na(pval) & pval < 0.05 & qval < 0.25)
  
  if(nrow(all_edges) == 0) {
    return(NULL)
  }
  
  all_edges <- all_edges %>% mutate(linetype = ifelse(coef > 0, "positive", "negative"), linesize = abs(coef))
  
  node_type <- data.frame(node = target_species, type = "Species") %>%
    rbind(data.frame(node = target_mets, type = "Metabolite")) %>%
    rbind(data.frame(node = target_clinicals, type = "Clinical"))
  
  active_nodes <- unique(c(all_edges$node1, all_edges$node2))
  node_table <- node_type %>% filter(node %in% active_nodes)
  
  g <- graph_from_data_frame(all_edges, vertices = node_table, directed = FALSE)
  E(g)$weight <- abs(all_edges$coef)
  g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = "first")
  g <- delete.vertices(g, which(degree(g) == 0))
  
  # ==========================================

  graphml_name <- sprintf("Network_Group_%s_Targeted14_P05_Q25.graphml", group_name)
  write_graph(g, graphml_name, format = "graphml")
  
  write.csv(node_table, sprintf("Nodes_Group_%s_Targeted14.csv", group_name), row.names = FALSE)
  write.csv(all_edges, sprintf("Edges_Group_%s_Targeted14.csv", group_name), row.names = FALSE)
  
}

# =========================================================
export_network_for_group("S")
export_network_for_group("C")