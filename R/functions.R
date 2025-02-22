#' @title generate_random_ps_subset
#' @description Generates a random subset of phyloseq object based on specified columns.
#' @param ps phyloseq object
#' @param level.columns column names based on which to subset
#' @return subsetted phyloseq object
#' @import phyloseq
#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by summarise across all_of n
#' @importFrom tidyr unite
#' @import magrittr 
generate_random_ps_subset <- function(ps, level.columns=NULL){
    if (is.null(level.columns) | (length(level.columns) == 0)) {
        sample.data.df <- sample_data(ps) %>% data.frame(check.names = FALSE)
        sample.names <- sample(rownames(sample.data.df), (nsamples(ps) * 95) %/% 100)
    } else {
        groups <- sample_data(ps) %>% 
            as_tibble(rownames = "id") %>% 
            as_tibble() %>% group_by(across(all_of(level.columns))) %>% 
            summarise(count = n(), ids = list(id), .groups='keep') %>%
            unite("group", all_of(level.columns), sep = "___", remove = FALSE)
        
        sample.names <- unlist(lapply(groups$ids, function(ids.list){
            n.all.samples.group <- length(ids.list)
            n.samples.min.group <- (n.all.samples.group * 80) %/% 100
            n.samples.max.group <- (n.all.samples.group * 95) %/% 100
            n.samples.batch.group <- sample(unlist(n.samples.min.group:n.samples.max.group), 1)
            sample(ids.list, n.samples.batch.group)
        }))
    }
    return(prune_samples(sample.names, ps))
}

#' @title bootstrap
#' @description Calculates pariwise correlations by bootstrapping
#' @param ps phyloseq object
#' @param level.columns column names based on which to separate into groups and calculate correlations
#' @param rho.pos.thr threshold which is used for filtering out pairs with lower positive correlation, Default: 0.8
#' @param rho.neg.thr threshold which is used for filtering out pairs with higher negative correlation, if not set -rho.pos.thr is taken, Default: NULL
#' @param p.adj.thr threshold which is used for filtering out pairs with larger p adjusted value, Default: 0.01
#' @param n.bootstrap.batches number of times to bootstrap, Default: 50
#' @param corr.pair.presence the ratio of bootstrap batches where the correlation pair is included, Default: 1
#' @return dataframe of correlation pairs
#' @importFrom dplyr filter mutate rowwise select summarize group_by n
#' @importFrom stats median p.adjust var
#' @import magrittr 
bootstrap <- function(ps, level.columns=NULL, rho.pos.thr=0.8, rho.neg.thr=NULL, p.adj.thr=0.01, n.bootstrap.batches=50, corr.pair.presence=1){
    ### takes ps, returns correlation table (does not do separately for categories)

    if (is.null(rho.neg.thr)){
        rho.neg.thr <- -rho.pos.thr
    } else if (rho.neg.thr > 0){
        rho.neg.thr <- -rho.neg.thr
    }

    if (!is.null(level.columns) && (length(level.columns) == 0)) {
        level.columns <- NULL
    }

    merged.corr.table <- data.frame()
    n.all.samples <- phyloseq::nsamples(ps)
    n.samples.min <- (n.all.samples * 85) %/% 100
    for (i in 1:n.bootstrap.batches){
        cat(paste0('iteration #', i, '\n'))
        ps.random.subset <- generate_random_ps_subset(ps, level.columns)
        corr.table <- phylosmith::co_occurrence(ps.random.subset,
                        treatment = level.columns) %>%
                        dplyr::filter((rho <= rho.neg.thr) | (rho >= rho.pos.thr)) %>%
                        mutate(p_adj = p.adjust(p, method = 'fdr')) %>%
                        filter(p_adj < p.adj.thr) %>%
                        rowwise() %>%
                        mutate(taxa_X = min(X, Y),
                                taxa_Y = max(X, Y)) %>%
                        select(-X, -Y) %>%
                        mutate(batch.id=i)

        merged.corr.table <- rbind(merged.corr.table, corr.table)
    }

    grouped.merged.corr.table <- merged.corr.table %>%
            group_by(taxa_X, taxa_Y, Treatment) %>%
            summarize(occurrence.batches.ratio = n() / n.bootstrap.batches,
                positive.rho = sum(rho > 0) / n(),
                negative.rho = sum(rho <= 0) / n(),
                p.adj = median(p_adj),
                rho = median(rho), 
                .groups = "keep"
            ) %>%
            filter(occurrence.batches.ratio > 0.8) %>%
            filter(positive.rho >= corr.pair.presence | negative.rho >= corr.pair.presence) %>%
            group_by(taxa_X, taxa_Y) %>%
            summarize(group.size.ratio = n() / length(unique(merged.corr.table$Treatment)),
                positive.rho = mean(positive.rho),
                negative.rho = mean(negative.rho),
                p.adj = median(p.adj),
                rho = median(rho),
                groups = paste(sort(Treatment), collapse = "___"),
                .groups = "keep"
            ) %>%
            filter(positive.rho == 1 | negative.rho == 1)
    return(grouped.merged.corr.table)
}


#' @title filter_ps
#' @description Filter the phyloseq object
#' @param ps phyloseq object
#' @param tax.lvl tax level, Default: 'Species'
#' @param relative.abundance make the abundance relative if set TRUE, Default: F
#' @param rabun.thr threshold which is used for filtering out abundances with lower relative abundance, Default: 0
#' @param prev.thr threshold which is used for filtering out abundances with lower prevalence, Default: 0
#' @param subset.columns named list for subsetting the phyloseq object based on sample data, where the names are column names and values are vectors of values to keep, Default: list()
#' @param level.columns columns based on which to separately the phyloseq object and separately do prevalence filtering, Default: NULL
#' @param species.names.modify if TRUE, changes the taxa_names names to consider the aggregation when tax.lvl is 'Species', Default: T
#' @return filtered phyloseq object
#' @import phyloseq
#' @importFrom dplyr filter mutate rowwise select summarize group_by inner_join n
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom stringr str_detect
#' @import magrittr
filter_ps <- function(ps, tax.lvl='Species', relative.abundance=F, rabun.thr=0, prev.thr=0, subset.columns=list(), level.columns=NULL, species.names.modify=T){
    if (!taxa_are_rows(otu_table(ps))) {
        otu_table_transposed <- otu_table(t(otu_table(ps)), taxa_are_rows = TRUE)
        ps@otu_table <- otu_table_transposed
    }

    filter.condition <- as.character(paste0(
        lapply(names(subset.columns), function(col) {
            vals <- subset.columns[[col]]
            paste0("(", col, " %in% c(", paste0(sprintf("'%s'", vals), collapse = ", "), "))")
        }),
        collapse = " & "
    ))

    if (filter.condition != ""){
        sample.data.ps <- data.frame(sample_data(ps))
        sample.data.ps.subset <- subset(sample.data.ps, eval(parse(text=filter.condition)))
        otu.table.ps.subset <- otu_table(ps)[, rownames(sample.data.ps.subset)]
        sample_data(ps) <- sample.data.ps.subset
        otu_table(ps) <- otu_table(as.matrix(otu.table.ps.subset), T)
        # ps <- subset_samples(ps, eval(parse(text = filter.condition)))
    }
    
    ps <- fantaxtic::name_na_taxa(ps, na_label = "<tax> <rank>")
    
    if (species.names.modify == T){
        aggregated.genera.table <- tax_table(ps) %>% data.frame() %>%
            group_by(Genus) %>% 
            summarise(
                size=n(),
                species.w.genus = sum(str_detect(Species, paste0(Genus, " Genus")))
            ) %>%
            filter((species.w.genus > 0) & (size > species.w.genus)) %>%
            mutate(
                species.w.genus.name = paste0(Genus, " Genus"),
                species.w.other.name = paste0(Genus, " Other")
            ) %>%
            tibble::column_to_rownames(var = "species.w.genus.name")
        
        species.names <- data.frame(tax_table(ps))[["Species"]]
        species.names <- ifelse(species.names %in% row.names(aggregated.genera.table), stringr::str_replace(species.names, 'Genus', 'Other'), gsub(" Genus", "", species.names))
        
        if (all(taxa_names(ps) == tax_table(ps)[, "Species"])){
            taxa_names(ps) <- species.names
        }
        
        tax_table(ps)[, "Species"]  <- species.names

        aggregated.genera.species <- tax_table(ps) %>% 
                            data.frame() %>% 
                            filter(Genus %in% aggregated.genera.table$Genus) %>% 
                            group_by(Genus) %>% 
                            summarise(Species = paste(Species, collapse = "___"))
    } else {
        aggregated.genera.species <- NULL
    }

    if (!is.null(tax.lvl)){
        ps <- ps %>%
            microViz::tax_fix() %>%
            microViz::tax_agg(tax.lvl) 
    }

    if (relative.abundance == T){
        cat(paste0('Filtering based on relative abundance ', rabun.thr, '\n'))
        ps <- transform_sample_counts(ps, function(x) x / sum(x))
        ps <- transform_sample_counts(ps, function(x) ifelse(x <= rabun.thr, 0, x))
    }
    
    if (prev.thr > 0) {
        cat(paste0('Filtering based on prevalence ', prev.thr, '\n'))
        sample.data.df <- sample_data(ps) %>% data.frame(check.names = FALSE) 
        if (!is.null(level.columns) && length(level.columns) > 0) {
            sample.data.df <- sample.data.df %>% select(level.columns)
        }
        otu.data <- otu_table(ps)
        if (taxa_are_rows(ps))  otu.data <- t(otu.data)
        otu.with.levels <- cbind(sample.data.df,  otu.data); otu.with.levels$id <- rownames(otu.with.levels)
        otu.with.levels.long <- otu.with.levels %>%
            pivot_longer(cols = colnames(otu.data), names_to = "Taxa", values_to = "Abundance") %>%
            mutate(Presence=Abundance>0)
        
        if (!is.null(level.columns) && length(level.columns) > 0) {
            prevalence.table <- otu.with.levels.long %>%
                group_by(across(all_of(level.columns)), Taxa) %>%
                summarize(Prevalence = mean(Presence), .groups = "drop") %>%
                filter(Prevalence > prev.thr)        
            otu.filtered <- otu.with.levels.long %>%
                inner_join(prevalence.table, by = c(level.columns, "Taxa")) %>%
                select(c('id', "Taxa", "Abundance")) %>%
                pivot_wider(names_from = "id", values_from = "Abundance", values_fill=0) %>%
                data.frame(check.names=FALSE)
        } else {
            prevalence.table <- otu.with.levels.long %>%
                group_by(Taxa) %>%
                summarize(Prevalence = mean(Presence), .groups = "drop") %>%
                filter(Prevalence > prev.thr)        
            otu.filtered <- otu.with.levels.long %>%
                inner_join(prevalence.table, by = c("Taxa")) %>%
                select(c('id', "Taxa", "Abundance")) %>%
                pivot_wider(names_from = "id", values_from = "Abundance", values_fill=0) %>%
                data.frame(check.names=FALSE)
        }
        row.names(otu.filtered) <- otu.filtered$Taxa; otu.filtered <- subset(otu.filtered, select=-c(Taxa))
        otu_table(ps) <- otu_table(otu.filtered, taxa_are_rows=T)
    }

    ps <- prune_taxa(taxa_sums(ps) > 0, ps)

    return(list(ps=ps, aggregated.genera.species=aggregated.genera.species))
}

#' @title assign_weights_lfc
#' @description Assigns discrete weights based on continuous values.
#' @param lfc.vec Vector of continuous values such as LFC. For discretization absolute values are conisdered.
#' @param quantiles Assigns weights based on quantiles, if set TRUE.
#' @param lfc.thr The vector of thresholds for discretization. Is treated as quantiles, if quantiles set TRUE, otherwise as absolute thresholds.
#' @param weights The vector of discrete weights for each possible case, Default: c(3, 1, -2)
#' @return The vector of corresponding weights.
#' @importFrom stats quantile
assign_weights_lfc <- function(lfc.vec, quantiles, lfc.thr, weights=c(3,1,-2)){
    if (length(lfc.thr) == length(weights) - 1 & length(lfc.thr) == 2){
        if (quantiles == T){
            thr.big <- quantile(abs(lfc.vec), lfc.thr[1], names=F, na.rm=T)
            thr.small <- quantile(abs(lfc.vec), lfc.thr[2], names=F, na.rm=T)
        }else{
            thr.big <- lfc.thr[1]
            thr.small <- lfc.thr[2] 
        }
        cat(paste0('The thresholds are ', round(thr.big, 3), ' ', round(thr.small, 3), '\n'))
        cat(paste0('The weights are ', weights[1], ' ', weights[2], ' ', weights[3], '\n'))
        lfc.vec.weighted <- sapply(lfc.vec, function(lfc){
            lfc.value <- abs(lfc)
            if (lfc.value >= thr.big){
                return(weights[1])
            }
            else if (lfc.value < thr.big & lfc.value >= thr.small){
                return(weights[2])
            }
            else{
                return(weights[3])
            }
        })
    } else if (length(lfc.thr) == length(weights) - 1 & length(lfc.thr) == 3){
        if (quantiles == T){
            thr.big <- quantile(abs(lfc.vec), lfc.thr[1], names=F)
            thr.mid <- quantile(abs(lfc.vec), lfc.thr[2], names=F)
            thr.small <- quantile(abs(lfc.vec), lfc.thr[3], names=F)
        }else{
            thr.big <- lfc.thr[1]
            thr.mid <- lfc.thr[2] 
            thr.small <- lfc.thr[3] 
        }
        cat(paste0('The thresholds are ', round(thr.big, 3), ' ', round(thr.mid, 2), ' ', round(thr.small, 3), '\n'))
        cat(paste0('The weights are ', weights[1], ' ', weights[2], ' ', weights[3], ' ', weights[4], '\n'))
        lfc.vec.weighted <- sapply(lfc.vec, function(lfc){
            lfc.value <- abs(lfc)

            if (lfc.value >= thr.big){
                return(weights[1])
            }
            else if (lfc.value < thr.big & lfc.value >= thr.mid){
                return(weights[2])
            }
            else if (lfc.value < thr.mid & lfc.value >= thr.small){
                return(weights[3])
            }else{
                return(weights[4])
            }
        })
    } else {
        cat('Provide appropriate thresholds.')
        return(NULL)
    }
    return(lfc.vec.weighted)
}


#' @title subgraph_solver_mwcsr
#' @description Returns the solved subgraph.
#' @param graph igraph object containing 'weights' attribute for vertices.
#' @return the solved subgraph
#' @import igraph 
subgraph_solver_mwcsr <- function(graph){
  rcsolver <- mwcsr::rmwcs_solver()
  graph <- igraph::upgrade_graph(graph)
  solved.subgraph <- mwcsr::solve_mwcsp(rcsolver, graph)
  return(solved.subgraph)
}


#' @title construct_network
#' @description Generates co-occurrence network.
#' @param ps phyloseq object
#' @param level.columns Vector of column names categories of which should be treated separately, for example c('Treatment'), Default: NULL
#' @param tax.lvl Taxonomic level, Default: 'Species'
#' @param subset.columns named list for subsetting the phyloseq object based on sample data, where the names are column names and values are vectors of values to keep, Default: list()
#' @param rabun.thr threshold which is used for filtering out abundances with lower relative abundance, Default: 0
#' @param prev.thr threshold which is used for filtering out abundances with lower prevalence, Default: 0
#' @param rho.pos.thr threshold which is used for filtering out pairs with lower positive correlation, Default: 0.8
#' @param rho.neg.thr threshold which is used for filtering out pairs with higher negative correlation, if not set -rho.pos.thr is taken, Default: NULL
#' @param p.adj.thr threshold which is used for filtering out pairs with larger p adjusted value, Default: 0.01
#' @param n.bootstrap.batches number of times to bootstrap, Default: 1000
#' @return named list, containing corr.table, network.object, aggregated.genera.species which are the correlation table, igraph oject and aggregated genera if tax.lvl is 'Species' respectively.
#' @export 
#' @importFrom dplyr mutate select
#' @importFrom stringr str_count
#' @importFrom stats setNames
#' @import magrittr
construct_network <- function(ps, level.columns=NULL, tax.lvl='Species', subset.columns=list(), rabun.thr=0, prev.thr=0, rho.pos.thr=0.8, rho.neg.thr=NULL, p.adj.thr=0.01, n.bootstrap.batches=1000){
    
    if (is.null(rho.neg.thr)){
        rho.neg.thr <- -rho.pos.thr
    } else if (rho.neg.thr > 0){
        rho.neg.thr <- -rho.neg.thr
    }

    filter.ps.output <- filter_ps(ps, tax.lvl=tax.lvl, relative.abundance=T, rabun.thr=rabun.thr, prev.thr=prev.thr, subset.columns=subset.columns, level.columns=level.columns)
    ps <- filter.ps.output$ps; aggregated.genera.species <- filter.ps.output$aggregated.genera.species

    grouped.merged.corr.table <- bootstrap(ps, level.columns, rho.pos.thr=rho.pos.thr, rho.neg.thr=rho.neg.thr, p.adj.thr=p.adj.thr, n.bootstrap.batches=n.bootstrap.batches)
    groups.categories <- unique(unlist(strsplit(grouped.merged.corr.table$groups, '___')))
    grouped.merged.corr.table <- grouped.merged.corr.table %>%
            mutate(across(
                all_of('groups'),
                list(!!!setNames(lapply(groups.categories, function(cat) { ~grepl(cat, ., ignore.case = TRUE)}), groups.categories)
                ), .names = "{fn}")
            ) %>%
            mutate(
                n.groups=str_count(groups, '___') + 1,
                group=ifelse(n.groups==1, groups, 'multiple')
            ) %>%
            mutate(absrho = abs(rho), rhosignpositive=(rho > 0)) %>%
            select(-c('group.size.ratio', 'positive.rho', 'negative.rho'))

    network.object <- igraph::graph_from_data_frame(grouped.merged.corr.table, directed=F)
    return(list(corr.table=grouped.merged.corr.table,
                network.object=network.object,
                aggregated.genera.species=aggregated.genera.species))
}


#' @title subgraph_analysis
#' @description Extract subgraph based on phyloseq object and co-occurrence network.
#' @param ps phyloseq object
#' @param network.object igraph object of the co-occurrence network
#' @param trt.column The treatment column name.
#' @param cntrl.level The category name for control in the trt.column.
#' @param trt.level The category name for treatment in the trt.column.
#' @param lfc.by.quantiles Assigns weights based on quantiles, if set TRUE.
#' @param lfc.thr The vector of thresholds for discretization. Is treated as quantiles, if quantiles set TRUE, otherwise as absolute thresholds.
#' @param lfc.weights The vector of discrete weights for each possible case.
#' @param tax.lvl Taxonomic level, Default: 'Species'
#' @param subset.columns Named list for subsetting the phyloseq object based on sample data, where the names are column names and values are vectors of values to keep, Default: list()
#' @return named list, containing igraph object of the network with subgraph attribute for vertices and edges, and ancombc results
#' @importFrom dplyr filter select rename mutate left_join sym
#' @importFrom igraph E V set_vertex_attr ends
#' @import phyloseq
#' @import magrittr 
subgraph_analysis <- function(ps, network.object, trt.column, cntrl.level, trt.level, lfc.by.quantiles, lfc.thr, lfc.weights, tax.lvl='Species', subset.columns=list()){
    cat(paste0('LFC thresholds: ', paste0(lfc.thr, collapse=' '), '\n'))
    cat(paste0('Weights: ', paste0(lfc.weights, collapse=' '),  '\n'))
    cat(paste0('LFC by quantiles: ', lfc.by.quantiles,  '\n'))
    subset.columns[[trt.column]] = c(trt.level, cntrl.level)

    ps <- fantaxtic::name_na_taxa(ps, na_label = "<tax> <rank>")
    ps.filtered <- filter_ps(ps, tax.lvl=NULL, relative.abundance=F, subset.columns=subset.columns)$ps
    ps.filtered <- phylosmith::set_treatment_levels(ps.filtered, treatment = trt.column, c(cntrl.level, trt.level))

    cat(paste0('Before filtering: nsamples = ', nsamples(ps), ', ntaxa = ', ntaxa(ps), '\n'))
    cat(paste0('After filtering: nsamples = ', nsamples(ps.filtered), ', ntaxa = ', ntaxa(ps.filtered), '\n'))

    lfc.column <- paste0("lfc_", trt.column, trt.level)
    network.vertices <- V(network.object)$name
    ancombc.raw <- ANCOMBC::ancombc2(verbose = F,
                        data = ps.filtered,
                        fix_formula = trt.column,
                        tax_level = tax.lvl,
                        p_adj_method = "fdr")$res 
    ancombc <- ancombc.raw %>%
                filter(!is.na(!!sym(lfc.column))) 
    
    cat(paste0('Number of lfc-s calculated: ', dim(ancombc)[1], '\n'))

    ancombc <- ancombc %>%
                filter(taxon %in% network.vertices)

    cat(paste0('Number of taxa with lfc-s present in the network: ', dim(ancombc)[1], '\n'))


    taxa.weighted <- ancombc %>%
        select(taxon, lfc.column) %>%
        rename(lfc = lfc.column) %>%
        mutate(weight = assign_weights_lfc(lfc, lfc.by.quantiles, lfc.thr, weights=lfc.weights))

    network.vertices.weighted <- data.frame(taxon = network.vertices) %>%
            left_join(taxa.weighted, by = "taxon") %>%
            mutate(weight = ifelse(is.na(weight), 0, weight))

    network.vertices.weighted <- network.vertices.weighted[match(V(network.object)$name, network.vertices.weighted$taxon),]

    weight.table <- table(network.vertices.weighted$weight)
    cat("Weight table:", paste(names(weight.table), weight.table, sep = ": ", collapse = ", "), '\n')
    
    network.object <- set_vertex_attr(network.object, 'weight', index = V(network.object), network.vertices.weighted$weight)
    solved.subgraph <- subgraph_solver_mwcsr(network.object)

    network.object <- set_vertex_attr(network.object, 'lfc', index = V(network.object), network.vertices.weighted$lfc)

    subgraph.vertices <- V(solved.subgraph$graph)$name
    vertex.attributes(network.object)$subgraph <- sapply(vertex.attributes(network.object)$name, function(x){
        return(x %in% subgraph.vertices)
    })

    E(network.object)$subgraph <- apply(ends(network.object, E(network.object)), 1, function(edge) {
        return((edge[1] %in% subgraph.vertices) & (edge[2] %in% subgraph.vertices))
    })

    return(list(network.object=network.object, ancombc=ancombc.raw))
}


#' @title plot_network
#' @description Plots co-occurrence network.
#' @param corr.table Dataframe of co-occurrence network, Default: NULL
#' @param network.object igraph object of the co-occurrence network, Default: NULL
#' @param taxa1.column.name The column name for the taxa 1 in the pairs, Default: 'taxa_X'
#' @param taxa2.column.name The column name for the taxa 2 in the pairs, Default: 'taxa_Y'
#' @param similarity.column.name The column name for similarity measure, Default: 'rho'
#' @importFrom igraph E V edge_attr graph_from_data_frame
#' @importFrom grDevices adjustcolor
#' @importFrom dplyr n
#' @importFrom graphics legend
#' @export
plot_network <- function(corr.table=NULL, network.object=NULL, taxa1.column.name='taxa_X', taxa2.column.name='taxa_Y', similarity.column.name='rho'){

    if (is.null(corr.table) & is.null(network.object)){
        print('Pass table or igraph network')
        return()
    }

    if (is.null(network.object)){
        network.object <- graph_from_data_frame(corr.table[c(taxa1.column.name, taxa2.column.name, similarity.column.name)], directed=F)
    }

    layout <- layout.fruchterman.reingold.grid

    
    plot(network.object, 
        layout = layout,
        vertex.label = V(network.object)$name, 
        vertex.size = 5, 
        edge.width = ifelse(edge_attr(network.object, similarity.column.name) > 0, 0.65, 0.85),
        edge.lty = ifelse(edge_attr(network.object, similarity.column.name)  > 0, 1, 2),
        vertex.color = V(network.object)$color,
        vertex.label.cex = 0.7, 
        vertex.label.color = adjustcolor("black", alpha.f = 0.9), 
        vertex.frame.color = adjustcolor("black", alpha.f = 0.2) 
    )


    legend("bottomleft", 
        legend = c("Positive Correlation", "Negative Correlation"), 
        col = c("black", "black"), 
        lty = c(1, 2), 
        lwd = 2, 
        cex = 0.59, 
        title = "Correlation Types"
    )
}


#' @title plot_subgraph
#' @description Plots subgraph network on co-occurrence network.
#' @param network.subgraph.object igraph object containing 'subgraph' attribute for vertices and edges.
#' @param similarity.column.name The column name for similarity measure, Default: 'rho'
#' @importFrom igraph layout.fruchterman.reingold.grid vertex.attributes V E
#' @importFrom grDevices adjustcolor
#' @importFrom tidyr replace_na
#' @importFrom graphics legend
#' @import magrittr 
#' @export
plot_subgraph <- function(network.subgraph.object, similarity.column.name='rho'){
    layout <- layout.fruchterman.reingold.grid
    lfc.colors.list = list("Positive LFC" = adjustcolor("darkolivegreen", alpha.f = 0.35), 
                            "Negative LFC" = adjustcolor("blue", alpha.f = 0.25), 
                            "NaN" = adjustcolor("grey", alpha.f = 0.4))
    vertex.attributes(network.subgraph.object)$color <- ifelse(
                            is.na(V(network.subgraph.object)$lfc),
                            lfc.colors.list[["NaN"]],
                            ifelse(
                                V(network.subgraph.object)$lfc > 0,
                                lfc.colors.list[["Positive LFC"]],
                                lfc.colors.list[["Negative LFC"]]
                            )
                        )

    plot(network.subgraph.object, 
        layout = layout,
        vertex.label = V(network.subgraph.object)$name, 
        vertex.size = 8 * abs(V(network.subgraph.object)$lfc %>% replace_na(min(abs(V(network.subgraph.object)$lfc), na.rm = TRUE))),
        edge.width = ifelse(E(network.subgraph.object)$subgraph, 1.1, 0.4),
        edge.color = ifelse(E(network.subgraph.object)$subgraph, "deeppink3", "grey"), 
        edge.lty = ifelse(E(network.subgraph.object)$rho > 0, 1, 2),
        vertex.color = V(network.subgraph.object)$color,
        vertex.label.cex = 0.35, 
        vertex.label.color = adjustcolor("black", alpha.f = 0.9),  
        vertex.frame.width = ifelse(V(network.subgraph.object)$subgraph, 2, 0.55),
        vertex.frame.color = ifelse(V(network.subgraph.object)$subgraph, adjustcolor("deeppink3", alpha.f = 0.9), adjustcolor("grey", alpha.f = 0.9))
    )

    legend("bottomright",
            legend = names(lfc.colors.list),
            col = unlist(unname(lfc.colors.list)),
            pch = 19,
            pt.cex = 1,
            cex = 0.59, 
            title = "Vertex categories")

    legend(x =-1.17, y = 1.1, 
        legend = c("Subgraph", "Other"),
        col = c("deeppink3", "grey"), 
        pch = 19,
        cex = 0.59,
        title = "Edge/Vertex Color")

    legend("bottomleft", 
        legend = c("Positive Correlation", "Negative Correlation"), 
        col = c("grey", "grey"), 
        lty = c(1, 2), 
        lwd = 2, 
        cex = 0.59, 
        title = "Correlation Types"
    )
}


#' @title rarefy
#' @description Rarefies phyloseq object.
#' @param ps phyloseq object
#' @param sample.sizes.list list of sample sizes with names of groups in level.column
#' @param level.column column name of the sample data based on which to rarefy the phyloseq object, i.e. Study
#' @return rarefied phyloseq object
#' @importFrom igraph layout.fruchterman.reingold.grid vertex.attributes V E
#' @import phyloseq
#' @importFrom microViz ps_filter
#' @importFrom dplyr pull sym
#' @import magrittr 
#' @export
rarefy <- function(ps, sample.sizes.list, level.column){
    categories <- sample_data(ps) %>%
                        pull(!!sym(level.column)) %>% 
                        unique() %>% 
                        unname()
    
    ps.subset.list <- c()
    for (category in categories){
        ps.subset <- ps_filter(ps, !!sym(level.column) == category)
        ps.subset <- rarefy_even_depth(ps.subset, sample.size=sample.sizes.list[[category]], rngseed=27091999)
        ps.subset.list <- c(ps.subset.list, ps.subset)
    }

    return(do.call(merge_phyloseq, ps.subset.list))
}


#' @title run.batches
#' @description run subgraph analysis for each batch
#' @param study study name
#' @param trt.column treatment column
#' @param cntrl.level control name
#' @param trt.level treatment name
#' @param experiment experiment name
#' @param ps.rarefied rarefied phyloseq object
#' @param network.object igraph network object of the co-occurrence network
#' @param tax.lvl Taxonomic level
#' @param lfc.thr The vector of thresholds for discretization. Is treated as quantiles, if quantiles set TRUE, otherwise as absolute thresholds.
#' @param lfc.weights The vector of discrete weights for each possible case.
#' @param lfc.by.quantiles Assigns weights based on quantiles, if set TRUE.
#' @return named list.
#' @importFrom igraph as_data_frame
#' @importFrom dplyr rename_with mutate
#' @import magrittr 
#' @export
run.batches <- function(study, trt.column, cntrl.level, trt.level, experiment, ps.rarefied, network.object, tax.lvl, lfc.thr, lfc.weights, lfc.by.quantiles){
    subset.columns <- list(Study=c(study)); subset.columns[[trt.column]] <- c(cntrl.level, trt.level)
    if (!is.null(experiment)) subset.columns[['Experiment']] <- c(experiment)
    suffix <- paste0(gsub(" ", "", study), '_', cntrl.level, '_', gsub(" ", "", trt.level))
    if (!is.null(experiment)) suffix <- paste0(suffix, '_', experiment)
    cat(paste0('Running ', suffix, '\n'))
    
    subgraph.analysis.output <- subgraph_analysis(ps.rarefied, network.object, trt.column=trt.column,
        cntrl.level=cntrl.level, trt.level=trt.level,
        tax.lvl=tax.lvl,
        lfc.by.quantiles=lfc.by.quantiles, lfc.thr=lfc.thr,
        lfc.weights=lfc.weights,
    #   lfc.by.quantiles=F, lfc.thr=c(2, 1),
        subset.columns=subset.columns)
    
    network.subgraph.object <- subgraph.analysis.output$network.object

    edge.atts <- colnames(as_data_frame(network.object, what = c("edges"))); vertex.atts <- colnames(as_data_frame(network.object, what = c("vertices")))

    subgraph.edges <- as_data_frame(network.subgraph.object, what = c("edges"))
    subgraph.edges.modified <- subgraph.edges %>%
                        rename_with(.fn = ~ paste0(.x, paste0('_', suffix)), .cols = -c(edge.atts))
    subgraph.vertices <- as_data_frame(network.subgraph.object, what = c("vertices")) %>%
                        mutate(abslfc = abs(lfc), lfcsignpositive=(lfc > 0)) 
    subgraph.vertices.modified <- subgraph.vertices %>%
                        rename_with(.fn = ~ paste0(.x, paste0('_', suffix)), .cols = -c(vertex.atts))
    return(list(subgraph.edges=subgraph.edges, subgraph.vertices=subgraph.vertices, 
                subgraph.edges.modified=subgraph.edges.modified,
                subgraph.vertices.modified=subgraph.vertices.modified,
                suffix=suffix, ancombc=subgraph.analysis.output$ancombc))
}


#' @importFrom stats dnorm var
.bias_em_patched <- function(beta, var_hat, tol, max_iter) {
  # Ensure beta and nu0 are filtered together
  neither_na = !(is.na(beta) | is.na(var_hat))
  beta = beta[neither_na]
  nu0 = var_hat[neither_na]
  
  if (any(nu0 == 0)) {
    stop_txt = sprintf(paste("Zero variances have been detected for the following taxa:", 
                             paste(names(which(nu0 == 0)), collapse = ", "), "Please remove these taxa or select a more parsimonious model", 
                             sep = "\n"))
    stop(stop_txt, call. = FALSE)
  }
  pi0_0 = 0.75
  pi1_0 = 0.125
  pi2_0 = 0.125
  delta_0 = mean(beta[beta >= quantile(beta, 0.25, na.rm = TRUE) & 
                        beta <= quantile(beta, 0.75, na.rm = TRUE)], na.rm = TRUE)
  if (is.na(delta_0)) 
    delta_0 = mean(beta, na.rm = TRUE)
  l1_0 = mean(beta[beta < quantile(beta, 0.125, na.rm = TRUE)], 
              na.rm = TRUE)
  if (is.na(l1_0)) 
    l1_0 = min(beta, na.rm = TRUE)
  l2_0 = mean(beta[beta > quantile(beta, 0.875, na.rm = TRUE)], 
              na.rm = TRUE)
  if (is.na(l2_0)) 
    l2_0 = max(beta, na.rm = TRUE)
  kappa1_0 = var(beta[beta < quantile(beta, 0.125, na.rm = TRUE)], 
                 na.rm = TRUE)
  if (is.na(kappa1_0) | kappa1_0 == 0) 
    kappa1_0 = 1
  kappa2_0 = var(beta[beta > quantile(beta, 0.875, na.rm = TRUE)], 
                 na.rm = TRUE)
  if (is.na(kappa2_0) | kappa2_0 == 0) 
    kappa2_0 = 1
  pi0_vec = pi0_0
  pi1_vec = pi1_0
  pi2_vec = pi2_0
  delta_vec = delta_0
  l1_vec = l1_0
  l2_vec = l2_0
  kappa1_vec = kappa1_0
  kappa2_vec = kappa2_0
  n_tax = length(beta)
  iterNum = 0
  epsilon = 100
  while (epsilon > tol & iterNum < max_iter) {
    pi0 = pi0_vec[length(pi0_vec)]
    pi1 = pi1_vec[length(pi1_vec)]
    pi2 = pi2_vec[length(pi2_vec)]
    delta = delta_vec[length(delta_vec)]
    l1 = l1_vec[length(l1_vec)]
    l2 = l2_vec[length(l2_vec)]
    kappa1 = kappa1_vec[length(kappa1_vec)]
    kappa2 = kappa2_vec[length(kappa2_vec)]
    pdf0 = vapply(seq(n_tax), function(i) dnorm(beta[i], 
                                                delta, sqrt(nu0[i])), FUN.VALUE = double(1))
    pdf1 = vapply(seq(n_tax), function(i) dnorm(beta[i], 
                                                delta + l1, sqrt(nu0[i] + kappa1)), FUN.VALUE = double(1))
    pdf2 = vapply(seq(n_tax), function(i) dnorm(beta[i], 
                                                delta + l2, sqrt(nu0[i] + kappa2)), FUN.VALUE = double(1))
    r0i = pi0 * pdf0/(pi0 * pdf0 + pi1 * pdf1 + pi2 * pdf2)
    r0i[is.na(r0i)] = 0
    r1i = pi1 * pdf1/(pi0 * pdf0 + pi1 * pdf1 + pi2 * pdf2)
    r1i[is.na(r1i)] = 0
    r2i = pi2 * pdf2/(pi0 * pdf0 + pi1 * pdf1 + pi2 * pdf2)
    r2i[is.na(r2i)] = 0
    pi0_new = mean(r0i, na.rm = TRUE)
    pi1_new = mean(r1i, na.rm = TRUE)
    pi2_new = mean(r2i, na.rm = TRUE)
    delta_new = sum(r0i * beta/nu0 + r1i * (beta - l1)/(nu0 + 
                                                          kappa1) + r2i * (beta - l2)/(nu0 + kappa2), na.rm = TRUE)/sum(r0i/nu0 + 
                                                                                                                          r1i/(nu0 + kappa1) + r2i/(nu0 + kappa2), na.rm = TRUE)
    l1_new = min(sum(r1i * (beta - delta)/(nu0 + kappa1), 
                     na.rm = TRUE)/sum(r1i/(nu0 + kappa1), na.rm = TRUE), 
                 0)
    if (is.na(l1_new)) 
      l1_new = 0
    l2_new = max(sum(r2i * (beta - delta)/(nu0 + kappa2), 
                     na.rm = TRUE)/sum(r2i/(nu0 + kappa2), na.rm = TRUE), 
                 0)
    if (is.na(l2_new)) 
      l2_new = 0
    obj_kappa1 = function(x) {
      log_pdf = log(vapply(seq(n_tax), function(i) dnorm(beta[i], 
                                                         delta + l1, sqrt(nu0[i] + x)), FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r1i * log_pdf, na.rm = TRUE)
    }
    kappa1_new = nloptr::neldermead(x0 = kappa1, fn = obj_kappa1, 
                                    lower = 0)$par
    obj_kappa2 = function(x) {
      log_pdf = log(vapply(seq(n_tax), function(i) dnorm(beta[i], 
                                                         delta + l2, sqrt(nu0[i] + x)), FUN.VALUE = double(1)))
      log_pdf[is.infinite(log_pdf)] = 0
      -sum(r2i * log_pdf, na.rm = TRUE)
    }
    kappa2_new = nloptr::neldermead(x0 = kappa2, fn = obj_kappa2, 
                                    lower = 0)$par
    pi0_vec = c(pi0_vec, pi0_new)
    pi1_vec = c(pi1_vec, pi1_new)
    pi2_vec = c(pi2_vec, pi2_new)
    delta_vec = c(delta_vec, delta_new)
    l1_vec = c(l1_vec, l1_new)
    l2_vec = c(l2_vec, l2_new)
    kappa1_vec = c(kappa1_vec, kappa1_new)
    kappa2_vec = c(kappa2_vec, kappa2_new)
    epsilon = sqrt((pi0_new - pi0)^2 + (pi1_new - pi1)^2 + 
                     (pi2_new - pi2)^2 + (delta_new - delta)^2 + (l1_new - 
                                                                    l1)^2 + (l2_new - l2)^2 + (kappa1_new - kappa1)^2 + 
                     (kappa2_new - kappa2)^2)
    iterNum = iterNum + 1
  }
  delta_em = delta_new
  pi1 = pi1_new
  pi2 = pi2_new
  l1 = l1_new
  l2 = l2_new
  kappa1 = kappa1_new
  kappa2 = kappa2_new
  C0 = which(beta >= quantile(beta, pi1, na.rm = TRUE) & beta < 
               quantile(beta, 1 - pi2, na.rm = TRUE))
  C1 = which(beta < quantile(beta, pi1, na.rm = TRUE))
  C2 = which(beta >= quantile(beta, 1 - pi2, na.rm = TRUE))
  nu = nu0
  nu[C1] = nu[C1] + kappa1
  nu[C2] = nu[C2] + kappa2
  wls_deno = sum(1/nu)
  wls_nume = 1/nu
  wls_nume[C0] = (wls_nume * beta)[C0]
  wls_nume[C1] = (wls_nume * (beta - l1))[C1]
  wls_nume[C2] = (wls_nume * (beta - l2))[C2]
  wls_nume = sum(wls_nume)
  delta_wls = wls_nume/wls_deno
  var_delta = 1/wls_deno
  if (is.na(var_delta)) 
    var_delta = 0
  output = c(delta_em = delta_em, delta_wls = delta_wls, var_delta = var_delta)
}


modified_co_occurrence <- function(
  phyloseq_obj,
  treatment    = NULL,
  subset       = NULL,
  rho          = 0,
  p            = 0.05,
  method       = "spearman",
  cores        = 1
) {
  phylosmith:::check_args(
    phyloseq_obj = phyloseq_obj,
    treatment    = treatment,
    subset       = subset,
    rho          = rho,
    p            = p,
    corr_method  = method,
    cores        = cores
  )

  if (cores == 0) cores <- (parallel::detectCores() - 1)
  phyloseq_obj <- 
    phylosmith::taxa_filter(phyloseq_obj, treatment, subset)
  treatment_name <- paste(treatment, collapse = phylosmith:::sep)

  treatment_classes <- 
    as.character(unique(phyloseq_obj@sam_data[[treatment_name]]))
  treatment_indices <- lapply(
    treatment_classes,
    FUN = function(trt) {
        which(as.character(phyloseq_obj@sam_data[[treatment_name]]) %in% trt)
    }
  )
  if (is.null(treatment)) {
    treatment_classes <- "Experiment_Wide"
    treatment_indices <- list(seq(nsamples(phyloseq_obj)))
  }
  too_few <- sapply(treatment_indices, length) < 3
  if (any(too_few)) {
    treatment_indices <- treatment_indices[!too_few]
  }
  phyloseq_obj <- phyloseq_obj@otu_table
  co_occurrence <- data.table::data.table()
  if (!is.vector(rho)) rho <- c(-rho, rho)
  for(i in seq_along(treatment_indices)){
    treatment_co_occurrence <- phylosmith:::Correlation(
      X               = phyloseq_obj[,treatment_indices[[i]]],
      cor_coef_cutoff        = rho,
      p_cutoff        = p,
      method          = method,
      ncores          = cores
    )
    treatment_co_occurrence[["X"]] <- rownames(
      phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[["X"]]]
    treatment_co_occurrence[["Y"]] <- rownames(
      phyloseq_obj[,treatment_indices[[i]]])[treatment_co_occurrence[["Y"]]]
    ### added
    if (dim(treatment_co_occurrence)[1] == 0){
      next
    }
    ###
    if(length(treatment_indices) > 0){
      treatment_co_occurrence <- cbind(
        Treatment = treatment_classes[i], 
        treatment_co_occurrence)
    }
    co_occurrence <- rbind(co_occurrence, treatment_co_occurrence)
  }
  
  return(data.table::as.data.table(co_occurrence))
}
