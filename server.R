server <- function(input, output) {
  
  codon_table <- data.frame(aa = as.data.frame(Biostrings::GENETIC_CODE)[,1], 
                            codon = Biostrings::GENETIC_CODE %>% as.data.frame() %>% rownames())
  codon_table[codon_table$aa == '*', 'aa'] <- 'Stop'
  deg_codons <-
    data.frame("M" = c("A", NA, NA, "C"),
               "V" = c("A", NA, 'G', "C"),
               "R" = c("A", NA, "G", NA),
               "H" = c("A", "T", NA, "C"),
               "W" = c("A", "T", NA, NA),
               "D" = c("A", "T", NA, "G"),
               "S" = c(NA, NA, "G", "C"),
               "B" = c(NA, "T", "G", "C"),
               "Y" = c(NA, "T", NA, "C"),
               "N" = c("A", "T", "G", "C"),
               "K" = c(NA, "T", "G", NA),
               "A" = c("A", NA, NA, NA),
               "T" = c(NA, "T", NA, NA),
               "G" = c(NA, NA, "G", NA),
               "C" = c(NA, NA, NA, "C")) %>% t() %>% data.frame()
  colnames(deg_codons) <- c("A", "T", "G", "C")
  
  if (file.exists('/home/ubuntu/AA2Cod/all_codons_df.csv')) {
    all_codons <- read_csv('/home/ubuntu/AA2Cod/all_codons_df.csv', show_col_types = FALSE) %>% as.data.frame()
  } else if (file.exists('all_codons_df.csv')) {
    all_codons <- read_csv('all_codons_df.csv', show_col_types = FALSE) %>% as.data.frame()
  }
  # } else {
  #   all_codons <- expand.grid(rownames(deg_codons), rownames(deg_codons), rownames(deg_codons))
  #   colnames(all_codons) <- c("pos1", 'pos2', 'pos3')
  #   all_codons <- cbind(all_codons, data.frame(matrix(nrow = nrow(all_codons), ncol = 21, dimnames = list(1:nrow(all_codons), c(Biostrings::sort(AA_STANDARD), '*')))))
  #   colnames(all_codons)[24] <- '*'
  #   for (row in 1:nrow(all_codons)) {
  #     cods <- expand.grid(as.character(deg_codons[all_codons[row,1],][!as.logical(is.na(deg_codons[all_codons[row,1],]))]), 
  #                         as.character(deg_codons[all_codons[row,2],][!as.logical(is.na(deg_codons[all_codons[row,2],]))]), 
  #                         as.character(deg_codons[all_codons[row,3],][!as.logical(is.na(deg_codons[all_codons[row,3],]))])
  #                         )
  #     cods$aa <- apply(cods,1, function(x) Biostrings::translate(DNAString(paste0(x, collapse = ''))) %>% as.character())
  #     countaas <- table(sort(cods$aa))
  #     all_codons[row, which(colnames(all_codons) %in% names(countaas))] <- countaas
  #   }
  #   write_csv(all_codons, 'all_codons_df.csv')
  # }
  all_codons_logic <- all_codons
  all_codons_logic[,] <- FALSE
  all_codons_logic[!is.na(all_codons)] <- TRUE
  all_codons_logic[,1:3] <- all_codons[,1:3]
  all_codons[is.na(all_codons)] <- 0
  all_aas <- apply(all_codons_logic,1, function(x) colnames(all_codons_logic)[4:24][as.logical(x[4:24])])
  names(all_aas) <- apply(all_codons_logic,1, function(x) paste0(x[1:3], collapse = ''))
  all_codons_logic$aas <- lapply(all_aas, function(x) paste0(x, collapse = '')) %>% as.character()
  all_codons$aas <- all_codons_logic$aas 
  
  
  data <- reactive({
    aa_selected <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "*")[c(input$A, input$C, input$D, input$E, input$F, input$G, input$H, input$I, input$K, input$L, input$M, input$N, input$P, input$Q, input$R, input$S, input$T, input$V, input$W, input$Y, input$Stop)]
    
    results <- data.frame(cods = names(all_aas), aas = all_codons$aas)
    results$length <- str_count(results$aas)
    results$intersect = lapply(all_aas, function(x) base::intersect(aa_selected, x) %>% length()) %>% unlist()
    results$extra <- results$length - results$intersect
    results$stop <- grepl('[*]', results$aas)
    
    results$total_cods <- apply(all_codons,1, function(x) x[str_split(x['aas'], '')[[1]]] %>% as.numeric() %>% sum())
    results$corr_cods <- apply(all_codons,1, function(x) x[aa_selected] %>% as.numeric() %>% sum())
    results$corr_cods_per <- round(results$corr_cods / results$total_cods * 100,2)
    results$wrong_cods <- results$total_cods - results$corr_cods
    
    results <- results[order(results$stop),]
    results <- results[order(results$extra),]
    results <- results[order(results$corr_cods_per, decreasing = TRUE),]
    results <- results[order(results$intersect, decreasing = TRUE),]
    if (length(aa_selected) > 0) {
      reg <- paste0('([^',paste0(aa_selected, collapse = ''),'])')
      results$aa <- str_replace_all(results$aas,reg,'<font color="red">\\1</font>')
    }
    
    cons_codon <- codon_table %>% filter(aa %in% aa_selected)
    cons_codon <- sapply(cons_codon$codon, function(x) str_split(x, '')) %>% as.data.frame() %>% t() %>% as.data.frame()
    pos1 <- apply(cons_codon, 2, unique)[1] %>% unlist()
    pos2 <- apply(cons_codon, 2, unique)[2] %>% unlist()
    pos3 <- apply(cons_codon, 2, unique)[3] %>% unlist()
    for (row in 1:nrow(deg_codons)) {
      if (length(pos1[!is.na(pos1)]) == length(deg_codons[row,][!is.na(deg_codons[row,])]) && all(sort(pos1) == sort(deg_codons[row,][!is.na(deg_codons[row,])]))) {
        deg1 <- rownames(deg_codons[row,])
      }
      if (length(pos2[!is.na(pos2)]) == length(deg_codons[row,][!is.na(deg_codons[row,])]) && all(sort(pos2) == sort(deg_codons[row,][!is.na(deg_codons[row,])]))) {
        deg2 <- rownames(deg_codons[row,])
      }
      if (length(pos3[!is.na(pos3)]) == length(deg_codons[row,][!is.na(deg_codons[row,])]) && all(sort(pos3) == sort(deg_codons[row,][!is.na(deg_codons[row,])]))) {
        deg3 <- rownames(deg_codons[row,])
      }
    }
    if (!exists("deg1")) {
      deg1 <- ''
    }
    if (!exists("deg2")) {
      deg2 <- ''
    }
    if (!exists("deg3")) {
      deg3 <- ''
    }
    
    codons_act <- expand.grid(deg_codons[deg1,][!is.na(deg_codons[deg1,])],
                              deg_codons[deg2,][!is.na(deg_codons[deg2,])],
                              deg_codons[deg3,][!is.na(deg_codons[deg3,])]) %>%
      apply(.,1, function(x) paste0(x, collapse = '')) %>% as.data.frame()
    codons_act$aa <- apply(codons_act, 1, function(x) codon_table[codon_table$codon == x['.'], 'aa'])
    colnames(codons_act) <- c('codon', 'aa')
    codon_des <- codon_table %>% filter(aa %in% aa_selected) %>% select(codon) %>% unlist()
    
    aa_des <- codon_table %>% filter(aa %in% aa_selected) %>% select(aa) %>% unlist() %>% unique() %>% sort()
    aa_act <- codons_act$aa %>% unique() %>% sort()
    
    list(aa_selected, pos1, pos2, pos3, deg1, deg2, deg3, codons_act, codon_des, aa_des, aa_act, results, p)
  })
  
  output$outcodontable <- renderDT({
    if (!identical(data()[[1]], character(0))) {
      out <- data()[[12]]
      out$aas <- NULL
      out[out$stop, 'stop'] <- "Yes"
      out[out$stop == FALSE, 'stop'] <- "No"
      colnames(out) <- c('Deg. Codon', '# of AAs', '# of correct AAs', '# of undesired AAs', '# Stop included', '# of codons', '# of correct codons', '% of correct codons', '# of wrong codons', 'amino acids')
      out
    }
  }, escape = FALSE, selection = 'single', server = FALSE, rownames = FALSE)
  
  output$aas <- renderPlot({
    if (!identical(data()[[1]], character(0)) && !is.null(input$outcodontable_rows_selected)) {
      aas <- str_split(gsub("<.*?>", "",data()[[12]][input$outcodontable_rows_selected, 'aa']),'')[[1]]
      cod <- str_split(data()[[12]][input$outcodontable_rows_selected, 'cods'], '')[[1]]
      p <- all_codons[all_codons$pos1 == cod[1] & all_codons$pos2 == cod[2] & all_codons$pos3 == cod[3],4:24]
      p <- pivot_longer(p, cols = colnames(p), names_to = 'AA', values_to = 'Count')
      p$AA <- factor(p$AA, levels = c("A","V","I","L","M","F","Y","W","D","E","K","R","H","S","T","N","Q","G","C","P", '*'))
      p$state <- 'Undesired'
      p[p$AA %in% data()[[1]], 'state'] <- 'Desired'
      ggplot(p) +
        geom_col(aes(x = AA, y = Count, fill = AA), color = 'black', size=0.3, width= 0.8, show.legend = F) +
        scale_y_continuous(expand = c(0,0), breaks = 1:10, limits = c(0, 6)) +
        facet_wrap(~state) +
        scale_fill_manual(
          values = c(
            "A" = 'grey90',
            "C" = '#82af8e',
            "D" = 'salmon',
            "E" = '#ff6969',
            "F"= 'grey40',
            "G"= 'white',
            "H"= '#119eff',
            "I"= 'grey52',
            "K"= '#6686ff',
            "L"= 'grey52',
            "M"= '#7b887e',
            "N"= '#ec93ea',
            "P"= '#94948a',
            "Q"= '#d1a3d0',
            "R"= 'blue',
            "S"= '#ff38f9',
            "T"= '#e86ae6',
            "V"= 'grey64',
            "W"= 'grey32',
            "Y"= '#6d5d6c',
            '*' = 'black'
          )
        ) +
        xlab('amino acid') +
        ylab('codon count') +
        theme_bw() +
        theme(axis.text = element_text(size=14), axis.title = element_text(size=15, face ='bold'), 
              panel.grid.minor.y = element_blank(), strip.text = element_text(size=15, face='bold'),
              panel.spacing = unit(2, "lines"))
    }
  })
  output$selectaa <- renderUI({
    if (all(!unlist(reactiveValuesToList(input)[c(AA_STANDARD, 'Stop')]))) {
      HTML('<i><font color="grey">Select an amino acid to see the a table for the best codons.</font></i>')
    } else {
      if (is.null(input$outcodontable_rows_selected)) {
        HTML('<i><font color="grey">Select a row in the data table to show codons.</font></i>')
      }
    }
  })
  
  observeEvent(input$allH, {
    lapply(c('A', 'F', 'W', 'I', 'L', 'Y', 'M', 'W', 'V'), function(x) updateMaterialSwitch(session = getDefaultReactiveDomain(), inputId = x, value = input$allH))
  })
  
  observeEvent(input$allC, {
    lapply(c('K', 'R', 'H', 'E', 'D'), function(x) updateMaterialSwitch(session = getDefaultReactiveDomain(), inputId = x, value = input$allC))
  })
  
  observeEvent(input$allS, {
    lapply(c('G', 'C', 'P', 'Stop'), function(x) updateMaterialSwitch(session = getDefaultReactiveDomain(), inputId = x, value = input$allS))
  })
  
  observeEvent(input$allP, {
    lapply(c('S', 'T', 'N', 'Q'), function(x) updateMaterialSwitch(session = getDefaultReactiveDomain(), inputId = x, value = input$allP))
  })
}