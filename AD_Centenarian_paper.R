#############################################################
# FINAL VERSION OF THE SCRIPT FOR THE NEW EFFECT-SIZE PAPER #
#############################################################

# LIBRARIES
    library(data.table)
    library(stringr)
    library(parallel)
    library(ggplot2)
    library(ggpubr)
    library(viridis)
    library(plyr)
    library(forestploter)
    library(grid)
    library(genpwr)
    library(nleqslv)

# FUNCTIONS
    # function to extract dosages and prepare for the prs
    function_extractAndPrepare <- function(ad_snps, data_path, pheno_path){
        all_dos = data.frame(bla = as.character())
        # grep snp id
        for (chr in unique(ad_snps$chrom)){
            snp_ids = c()
            tmp = ad_snps[which(ad_snps$chrom == chr), ]
            for (i in 1:nrow(tmp)){
                tmp_grep = system(paste0('grep -w ', tmp$pos[i], ' ', data_path, '/chr', chr, '.dose.unscrambled.pvar'), intern = T)
                tmp_grep = data.frame(str_split_fixed(tmp_grep, '\t', 7), stringsAsFactors=F)
                snp_ids = c(snp_ids, tmp_grep$X3)
            }
            # write the list of ids to grep
            write.table(snp_ids, paste0('chr', chr, '_snp_ids.txt'), quote = F, row.names = F, col.names = F)
            # then extract dosages
            cmd = paste0('plink2 --pfile ', data_path, '/chr', chr, '.dose.unscrambled --extract chr', chr, '_snp_ids.txt --export A --out chr', chr, '_dosages')
            system(cmd)
            # read dosages
            tmp_dos = fread(paste0('chr', chr, '_dosages.raw'), h=T, stringsAsFactors=F)
            x <- colnames(tmp_dos); col.index <- grep(":", x); dos_clean <- tmp_dos[, ..col.index]; dos_clean$IID <- tmp_dos$IID
            # sort by id
            dos_clean = dos_clean[order(dos_clean$IID), ]
            if (nrow(all_dos) == 0){
                all_dos = dos_clean
            } else {
                all_dos = merge(all_dos, dos_clean, by = 'IID')
            }
        }

        # add pheno
        dosages_pheno <- merge(all_dos, pheno_path, by.x="IID", by.y = 'ID_GWAS')

        return(dosages_pheno)
    }

    # function to calculate PRS given the full list of snp dosages, the snp information from GWAS, and the name of the phenotype to take the corresponding SNPs
    function_PRS <- function(dosages_pheno, ad_snps){
        snplist = colnames(dosages_pheno)[grep('chr', colnames(dosages_pheno))]
        cat(paste0('** unique snps included for this phenotype --> ', length(snplist), '\n'))

        # make prs
        prs <- as.data.frame(matrix(data=0, nrow=nrow(dosages_pheno), ncol=3))
        colnames(prs) <- c("IID", "PRS", "PHENO")
        prs$IID <- dosages_pheno$IID
        prs$PHENO <- dosages_pheno$diagnosis
        prs$PC1 = dosages_pheno$PC1
        prs$PC2 = dosages_pheno$PC2
        prs$PC3 = dosages_pheno$PC3
        prs$PC4 = dosages_pheno$PC4
        prs$PC5 = dosages_pheno$PC5

        # main loop
        for (snp in snplist){
            # check if we need to include the snp
            snp_name = snp
            snp_pos = str_split_fixed(snp_name, ':', 4)[, 2]
            snp_allele = str_split_fixed(snp_name, '_', 2)[, 2]
            if (snp_pos %in% ad_snps$pos){
                # select snp
                snp_data <- as.vector(unlist(dosages_pheno[, ..snp]))
                # extract snp information
                tmp_snp_info = ad_snps[which(ad_snps$pos == snp_pos), ]
                effect_allele_gwas = str_split_fixed(tmp_snp_info$"minor/major", "/", 2)[, 1]
                beta_gwas = log(as.numeric(tmp_snp_info$or))
                # check allele and in case flip dosages
                if (snp_allele != effect_allele_gwas){ snp_data = 2 - snp_data }
                # get tmp score
                tmp <- beta_gwas * snp_data
                prs$PRS <- prs$PRS + tmp
            }
        }
        return(prs)
    }

    # function to test for association
    function_testAssoc = function(prs, pheno, pheno_final_raw, do_pca){
        # first, scale the prs
        prs$PRS = scale(prs$PRS)
        # isolate the children
        children = prs[which(prs$PHENO %in% c('family_100plus')),]
        # keep only samples that passed qc and add pcs
        pcs = pheno[, c('ID_GWAS', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
        prs = merge(prs, pcs, by.x = 'IID', by.y = 'ID_GWAS')
        # chc vs. controls
        tmp_prs = prs[which(prs$PHENO %in% c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO == 'Centenarian')] = 0; print(table(tmp_prs$STATUS))
        model_chc_ctr = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_chc_ctr = as.data.frame(t(summary(model_chc_ctr)$coefficients[2, ]))
        res_chc_ctr$comparison = 'chc_ctr'
        res_chc_ctr$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_chc_ctr$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # chc vs. ad
        tmp_prs = prs[which(prs$PHENO %in% c("Centenarian", "Probable_AD", "Possible_AD", "AD_path")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO == 'Centenarian')] = 0; print(table(tmp_prs$STATUS))
        model_chc_ad = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_chc_ad = as.data.frame(t(summary(model_chc_ad)$coefficients[2, ]))
        res_chc_ad$comparison = 'chc_ad'
        res_chc_ad$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_chc_ad$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # ad vs. ctr
        tmp_prs = prs[which(prs$PHENO %in% c("Control_LASA", "Control_other_twin", "Control_100plus", "Control_path", "SCD", "Probable_AD", "Possible_AD", "AD_path")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 0; print(table(tmp_prs$STATUS))
        model_ad_ctr = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_ad_ctr = as.data.frame(t(summary(model_ad_ctr)$coefficients[2, ]))
        res_ad_ctr$comparison = 'ad_ctr'
        res_ad_ctr$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_ad_ctr$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # finally children vs. their partners
        children$STATUS = 0
        # what if we compare children with the controls
        tmp_prs = rbind(children[which(children$STATUS == 0), c('IID', 'PRS', 'PHENO')], prs[which(prs$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")), c('IID', 'PRS', 'PHENO')])
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO == 'family_100plus')] = 0
        # we need to check the childrens
        if (do_pca == 'no_pca'){
            pcs = fread('pca_children_controls.eigenvec', h=F, stringsAsFactors=F)
            colnames(pcs) = c('FID', 'IID', paste0('PC', seq(1, 20)))
            pcs_sb = pcs[, c('IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
            children_clean = merge(tmp_prs, pcs_sb, by = 'IID')
        } else {
            set.seed(020202)
            children_clean = checkchildren(tmp_prs, pheno_final_raw)
        }
        # model
        model_children = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = children_clean, family = 'binomial')
        summary(model_children)
        res_children = as.data.frame(t(summary(model_children)$coefficients[2, ]))
        res_children$comparison = 'children_vs_ctr'
        res_children$N_cases = nrow(children_clean[which(children_clean$STATUS == 1),]); res_children$N_ctr = nrow(children_clean[which(children_clean$STATUS == 0),])
        # put results together
        results = rbind(res_chc_ad, res_ad_ctr, res_chc_ctr, res_children)
        return(results)
    }

    # generic function to run single-variant association tests
    function_singleVar <- function(dosages_pheno, pheno){
        snp_assoc = list()
        snplist = colnames(dosages_pheno)[grep(':', colnames(dosages_pheno))]
        for (snp in 1:length(snplist)){
            cat(paste0('** processing snp ', snp, '/', length(snplist), '    \r'))
            snp = snp + 1
            # check if we need to include the snp
            snp_name = colnames(dosages_pheno)[snp]
            snp_allele = str_split_fixed(snp_name, '_', 2)[, 2]
            # select snp
            snp_data <- data.frame(iid = dosages_pheno$IID, dos = as.vector(unlist(dosages_pheno[, ..snp])), phenotype = dosages_pheno$diagnosis, stringsAsFactors=F)
            # set binary phenotypes for chc vs. controls
            snp_data_1 = snp_data[which(snp_data$phenotype %in% c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
            snp_data_1$phenotype = as.character(snp_data_1$phenotype)
            snp_data_1$phenotype_binary = 1; snp_data_1$phenotype_binary[which(snp_data_1$phenotype == 'Centenarian')] = 0
            snp_data_1 = merge(snp_data_1, pheno, by.x = 'iid', by.y = 'ID_GWAS')
            model_1 = glm(phenotype_binary ~ dos + PC1 + PC2 + PC3 + PC4 + PC5, data = snp_data_1, family = 'binomial')
            df_chc_ctr = data.frame(t(summary(model_1)$coefficients[2,])); df_chc_ctr$snp_name = snp_name
            colnames(df_chc_ctr) = c('beta', 'se', 'z_value', 'p', 'snp')
            df_chc_ctr$test = 'ctr_vs_chc'
            # set binary phenotypes for chc vs. AD
            snp_data_1 = snp_data[which(snp_data$phenotype %in% c("Centenarian", "Probable_AD", "Possible_AD", "AD_path")),]
            snp_data_1$phenotype = as.character(snp_data_1$phenotype)
            snp_data_1$phenotype_binary = 1; snp_data_1$phenotype_binary[which(snp_data_1$phenotype == 'Centenarian')] = 0
            snp_data_1 = merge(snp_data_1, pheno, by.x = 'iid', by.y = 'ID_GWAS')
            model_1 = glm(phenotype_binary ~ dos + PC1 + PC2 + PC3 + PC4 + PC5, data = snp_data_1, family = 'binomial')
            df_chc_ad = data.frame(t(summary(model_1)$coefficients[2,])); df_chc_ad$snp_name = snp_name
            colnames(df_chc_ad) = c('beta', 'se', 'z_value', 'p', 'snp')
            df_chc_ad$test = 'ad_vs_chc'
            # set binary phenotypes for chc vs. AD
            snp_data_1 = snp_data[which(snp_data$phenotype %in% c("Probable_AD", "Possible_AD", "AD_path", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
            snp_data_1$phenotype = as.character(snp_data_1$phenotype)
            snp_data_1$phenotype_binary = 0; snp_data_1$phenotype_binary[which(snp_data_1$phenotype %in% c("Probable_AD", "Possible_AD", "AD_path"))] = 1
            snp_data_1 = merge(snp_data_1, pheno, by.x = 'iid', by.y = 'ID_GWAS')
            model_1 = glm(phenotype_binary ~ dos + PC1 + PC2 + PC3 + PC4 + PC5, data = snp_data_1, family = 'binomial')
            df_ctr_ad = data.frame(t(summary(model_1)$coefficients[2,])); df_ctr_ad$snp_name = snp_name
            colnames(df_ctr_ad) = c('beta', 'se', 'z_value', 'p', 'snp')
            df_ctr_ad$test = 'ad_vs_ctr'
            # combine associations
            df = rbind(df_chc_ctr, df_chc_ad, df_ctr_ad)
            df$allele = snp_allele
            snp_assoc[[(length(snp_assoc) + 1)]] = df
        }
        snp_assoc = rbindlist(snp_assoc)
        return(snp_assoc)
    }

    # effect size ratio
    effect_size_ratio <- function(singleAssoc, ad_snps){
        # output
        ratios = list()
        # derive list of unique snps
        snps = unique(singleAssoc$snp)
        # add locus to ad_snps
        ad_snps$locus = paste(ad_snps$chrom, ad_snps$pos, sep = ":")
        # also add standard error
        tmp = data.frame(str_split_fixed(ad_snps$ci, '–|-', 2), stringsAsFactors = F)
        ad_snps$se_gwas = abs((log(as.numeric(tmp$X1)) - log(as.numeric(ad_snps$or))) / 1.96)
        # loop on snps
        for (s in snps){
            # get snp locus (chromosome:position)
            snp_locus = paste(str_split_fixed(s, ':', 4)[, 1], str_split_fixed(s, ':', 4)[, 2], sep = ':'); snp_locus = str_replace_all(snp_locus, 'chr', '')
            # take gwas associations
            gwas_beta = log(as.numeric(ad_snps$or[which(ad_snps$locus == snp_locus)]))
            gwas_allele = str_split_fixed(ad_snps$"minor/major"[which(ad_snps$locus == snp_locus)], '/', 2)[, 1]
            # gather my associations
            my_assoc = singleAssoc[which(singleAssoc$snp == s),]
            my_allele = unique(my_assoc$allele)
            # check alleles -- align with gwas --> in case, flip my beta (beta_aligned)
            my_assoc$gwas_allele = gwas_allele; my_assoc$gwas_beta = gwas_beta; my_assoc$beta_aligned = my_assoc$beta
            if (gwas_allele != my_allele){ my_assoc$beta_aligned = my_assoc$beta_aligned * (-1) }
            # save
            ratios[[(length(ratios) + 1)]] = my_assoc
        }
        ratios = rbindlist(ratios)
        # no need to standardize now!
        # finally calculate ratios
        ratios$ratio_beta = ratios$beta_aligned / ratios$gwas_beta
        # add gene information
        ratios$locus = paste(str_split_fixed(ratios$snp, ':', 4)[, 1], str_split_fixed(ratios$snp, ':', 4)[, 2], sep = ':'); ratios$locus = str_replace_all(ratios$locus, 'chr', '')
        gene_info = ad_snps[, c('locus', 'gene', 'se_gwas')]; ratios = merge(ratios, gene_info, by = 'locus')
        # add the sampling framework
        ratios$one_tail_p = NA; ratios$two_tail_p = NA
        ratios$low_ci_ratio = NA; ratios$up_ci_ratio = NA
        for (i in 1:nrow(ratios)){
            # non-standardized
            tmp_my <- rnorm(n=10000, mean = ratios$beta_aligned[i], sd = ratios$se[i])
            tmp_gw <- rnorm(n=10000, mean = ratios$gwas_beta[i], sd = ratios$se_gwas[i])
            tmp_ratios <- tmp_my / tmp_gw
            # calculate up and low ci for standardized and non-standardized
            tmp_low = quantile(tmp_ratios, .025); tmp_up = quantile(tmp_ratios, .975)
            cdf <- ecdf(tmp_ratios); cdf2 <- cdf(1)
            if (cdf2 > 0.5){ ratios$one_tail_p[i] = 1 - cdf2; ratios$two_tail_p[i] = ratios$one_tail_p[i] * 2 } else { ratios$one_tail_p[i] = cdf2; ratios$two_tail_p[i] = ratios$one_tail_p[i] * 2 }
            ratios$low_ci_ratio[i] = tmp_low; ratios$up_ci_ratio[i] = tmp_up
        }

        return(ratios)
    }
    
    # function to check significance, correct for multiple tests, check directions and binomial tests -- v2
    checkAssociations <- function(ratios_tmp, ad_snps){
        # correct for multiple tests using fdr
        ratios_tmp$p_adjust = p.adjust(ratios_tmp$p, 'fdr')
        # show number of nominal significant hits
        cat(paste0('** # of nominal significant hits --> ', nrow(ratios_tmp[which(ratios_tmp$p <= 0.05),]), '\n'))
        # show number of significant hits after correction
        cat(paste0('** # of significant hits after fdr correction --> ', nrow(ratios_tmp[which(ratios_tmp$p_adjust <= 0.05),]), '\n'))
        # check direction
        ratios_tmp$direction = ifelse((ratios_tmp$gwas_beta * ratios_tmp$beta_aligned) >0, 'same', 'different')
        cat(paste0('** ', nrow(ratios_tmp[which(ratios_tmp$direction == 'same'),]), '/', nrow(ratios_tmp), ' SNPs are in the same direction\n'))
        print(binom.test(x = nrow(ratios_tmp[which(ratios_tmp$direction == 'same'),]), n = nrow(ratios_tmp), p = 0.5))
        # average change in effect size
        cat(paste0('** Avg. change in effect size (non-standardized) is ', round(median(ratios_tmp$ratio_beta), 2), ' with a deviation of ', round(sd(ratios_tmp$ratio_beta), 2), '\n'))
        # number of snps with change >1
        cat(paste0('** # of SNPs with ratio >1 is --> ', nrow(ratios_tmp[which(ratios_tmp$ratio_beta >= 1),])), '\n')
        cat(paste0('** # of SNPs with ratio >0 is --> ', nrow(ratios_tmp[which(ratios_tmp$ratio_beta >= 0),])), '\n')
        cat('\n\n')
    }

    # function to extract frequencies for the maf plot
    extract_frequencies <- function(ratios){
        snps_interest = str_split_fixed(unique(ratios$snp), '_', 2)[, 1]
        tmp = str_split_fixed(snps_interest, ':', 4)
        df = data.frame(snpid = snps_interest, chrom = tmp[, 1], position = tmp[, 2])
        # also need to write samples to look at
        pheno <- read.table("/project/holstegelab/Share/gwas_array/mapping_files/20211027_phenotypes_ADC_samplesKept.txt", h=T)
        chc = pheno[which(pheno$diagnosis == 'Centenarian'),]
        ctr = pheno[which(pheno$diagnosis %in% c("Control_100plus", "Control_LASA", "Control_other", "Control_other_twin", "Control_path", "SCD")),]
        ad = pheno[which(pheno$diagnosis %in% c("Probable_AD", "Possible_AD", "AD_path")),]
        write.table(chc$ID_GWAS, 'chc_samples.txt', quote=F, row.names = F, col.names = F)
        write.table(ctr$ID_GWAS, 'ctr_samples.txt', quote=F, row.names = F, col.names = F)
        write.table(ad$ID_GWAS, 'ad_samples.txt', quote=F, row.names = F, col.names = F)
        # main loop to calculate frequencies
        all_frequencies = list()
        all_frequencies_plot = list()
        for (chr in unique(df$chrom)){
            tmp_snps = df[which(df$chrom == chr),]
            write.table(tmp_snps$snpid, paste0(chr, '_snps.txt'), quote=F, row.names=F, col.names=F)
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep chc_samples.txt --extract ', chr, '_snps.txt --freq --out ', chr, '_freq_chc'))
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep ad_samples.txt --extract ', chr, '_snps.txt --freq --out ', chr, '_freq_ad'))
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep ctr_samples.txt --extract ', chr, '_snps.txt --freq --out ', chr, '_freq_ctr'))
            # read frequencies back in
            chc_freq = fread(paste0(chr, '_freq_chc.afreq'), h=T)
            ctr_freq = fread(paste0(chr, '_freq_ctr.afreq'), h=T)
            ad_freq = fread(paste0(chr, '_freq_ad.afreq'), h=T)
            # combine frequencies
            df_freq = data.frame(snpid = chc_freq$ID, ref = chc_freq$REF, alt = chc_freq$ALT, n_allele_chc = chc_freq$OBS_CT, n_allele_ctr = ctr_freq$OBS_CT, n_allele_ad = ad_freq$OBS_CT, freq_chc = chc_freq$ALT_FREQS, freq_ctr = ctr_freq$ALT_FREQS, freq_ad = ad_freq$ALT_FREQS)
            df_freq_plot = rbind(chc_freq, ctr_freq, ad_freq); df_freq_plot$Phenotype = c(rep('CHC', nrow(chc_freq)), rep('CTR', nrow(ctr_freq)), rep('AD', nrow(ad_freq)))
            # save
            all_frequencies[[(length(all_frequencies) + 1)]] = df_freq
            all_frequencies_plot[[(length(all_frequencies_plot) + 1)]] = df_freq_plot
        }
        all_frequencies = rbindlist(all_frequencies)
        tmp = str_split_fixed(all_frequencies$snpid, ':', 4)
        all_frequencies$locus = paste0(tmp[, 1], ':', tmp[, 2]); all_frequencies$locus = str_replace_all(all_frequencies$locus, 'chr', '')
        all_frequencies_plot = rbindlist(all_frequencies_plot)
        tmp2 = str_split_fixed(all_frequencies_plot$ID, ':', 4)
        all_frequencies_plot$locus = paste0(tmp2[, 1], ':', tmp2[, 2]); all_frequencies_plot$locus = str_replace_all(all_frequencies_plot$locus, 'chr', '')
        res = list(all_frequencies, all_frequencies_plot)
        return(res)
    }

    # function to draw figure 1 -- 3 plots: maf, effect-size and ratio
    plot_figure1 <- function(singleAssoc, ratios, ad_snps, all_freqs_plot, snps_info){
        # create right order for the snps
            tmp = ratios[which(ratios$test == 'ad_vs_chc')]; tmp = tmp[order(tmp$ratio_beta),]; 
            tmp$gene[which(duplicated(tmp$gene) == TRUE)] = paste0(tmp$gene[which(duplicated(tmp$gene) == TRUE)], ' (2)')
            tmp$gene = factor(tmp$gene, levels = tmp$gene)
            ad_snps_new_info = ad_snps[, c('locus', 'new_known')]
            tmp = merge(tmp, ad_snps_new_info, by = 'locus')
        # plot 1 is the MAF plot
            # finally let's do a plot
            all_freqs_plot = merge(all_freqs_plot, snps_info, by = 'locus')
            all_freqs_plot$MAF_FREQS = ifelse(all_freqs_plot$ALT_FREQS < 0.5, all_freqs_plot$ALT_FREQS, 1 - all_freqs_plot$ALT_FREQS)
            # reorder before the plot
            ratios_values = tmp[, c('ratio_beta', 'locus', 'gene', 'new_known')]; all_freqs_plot = merge(all_freqs_plot, ratios_values, by = 'locus')
            all_freqs_plot = all_freqs_plot[order(all_freqs_plot$ratio_beta),]
            # change labels and we're done
            all_freqs_plot$Phenotype[which(all_freqs_plot$Phenotype == 'AD')] = 'AD cases'
            all_freqs_plot$Phenotype[which(all_freqs_plot$Phenotype == 'CTR')] = 'Healthy controls'
            all_freqs_plot$Phenotype[which(all_freqs_plot$Phenotype == 'CHC')] = 'Centenarians'
            # adjust labels for new and known snps
            all_freqs_plot$label = ifelse(all_freqs_plot$new_known == 'new', 'New SNP in latest GWAS', 'Known SNP')
            # plot
            plt = ggplot(all_freqs_plot, aes(x = gene, y = MAF_FREQS, color = Phenotype, shape = Phenotype)) + geom_point(size = 3) + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + 
                xlab('') + ylab('Minor Allele Frequency')
            plt = plt + facet_grid(cols = vars(label), drop = FALSE, scales = "free")
        # plot 2 is the grouped barplot with effect-sizes
            # prepare data for grouped barplot
            grp_plot = list()
            for (i in 1:nrow(tmp)){ grp_plot[[(length(grp_plot) + 1)]] = data.frame(Gene = rep(tmp$gene[i], 2), Effect = c(tmp$gwas_beta[i], tmp$beta_aligned[i]), Study = c('GWAS', 'AD cases vs. Centenarians'), low = c(tmp$gwas_beta[i] - (1.96*tmp$se_gwas[i]), tmp$beta_aligned[i] - (1.96*tmp$se[i])), up = c(tmp$gwas_beta[i] + (1.96*tmp$se_gwas[i]), tmp$beta_aligned[i] + (1.96*tmp$se[i]))) }
            grp_plot = rbindlist(grp_plot)
            grp_plot[which(grp_plot$Gene == 'TREML2' & grp_plot$Study != 'GWAS'), c('Effect', 'low', 'up')] = NA
            # add labels for new of known snp
            new_snps = all_freqs_plot$gene[which(all_freqs_plot$new_known == 'new')]
            grp_plot$label = ifelse(grp_plot$Gene %in% new_snps, 'New SNP in latest GWAS', 'Known SNP')
            # plot
            plt_grp = ggplot(grp_plot, aes(fill=Study, y=Effect, x=Gene)) + geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim = c(-1, 1), clip = "on") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + geom_errorbar(aes(ymin = low, ymax = up), width = 0.2, color = 'grey60', position=position_dodge(.9)) + ylim(-1, 1)
            plt_grp = plt_grp + facet_grid(cols = vars(label), drop = FALSE, scales = "free")
        # plot 3 is the effect-size ratio
            tmp$Comparison_symbol = ''; tmp$Comparison_symbol[which(tmp$two_tail_p <= 0.05)] = 'x'
            tmp$Comparison = 'Not different'; tmp$Comparison[which(tmp$two_tail_p <= 0.05)] = 'Different from published'
            tmp$adj_p = p.adjust(tmp$p, 'fdr')
            # add labels for new of known snp
            tmp$label = ifelse(tmp$gene %in% new_snps, 'New SNP in latest GWAS', 'Known SNP')
            # find significant associations to be annotated
            star_df = data.frame(x = tmp$gene[which(tmp$p <= 0.05)], y = tmp$up_ci_ratio[which(tmp$p <= 0.05)], comparison = tmp$Comparison[which(tmp$p <= 0.05)], label = tmp$label[which(tmp$p <= 0.05)])
            fdr_df = data.frame(x = tmp$gene[which(tmp$adj_p <= 0.05)], y = tmp$up_ci_ratio[which(tmp$adj_p <= 0.05)], comparison = tmp$Comparison[which(tmp$adj_p <= 0.05)], label = tmp$label[which(tmp$adj_p <= 0.05)])
            tmp$ratio_beta[which(tmp$gene == 'TREML2')] = NA; tmp$low_ci_ratio[which(tmp$gene == 'TREML2')] = NA; tmp$up_ci_ratio[which(tmp$gene == 'TREML2')] = NA
            # the actual plot
            plt_ratio = ggplot(tmp, aes(x = gene, y = ratio_beta, fill = Comparison)) + geom_bar(stat = 'identity') + geom_errorbar(aes(ymin = low_ci_ratio, ymax = up_ci_ratio), width = 0.4, color = 'grey60') + geom_hline(yintercept = 1, linetype = 'dashed', col = 'red') +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + xlab('') + ylab('Effect-size Ratio')
            plt_ratio = plt_ratio + facet_grid(cols = vars(label), drop = FALSE, scales = "free")
            # annotation of NA (TREML2)
            ann_text_treml2 = data.frame(gene = 'TREML2', ratio_beta = 0, label = 'Known SNP', Comparison = 'Different from published')
            plt_ratio = plt_ratio + geom_text(data = ann_text_treml2, label = 'NA', size = 2)
            # annotation of significant associations
            ann_sign_assoc = data.frame(gene = star_df$x, ratio_beta = star_df$y, Comparison = star_df$comparison, label = star_df$label)
            plt_ratio = plt_ratio + geom_text(data = ann_sign_assoc, label = '*')
            # annotation of fdr significant associations
            ann_fdr_assoc = data.frame(gene = fdr_df$x, ratio_beta = fdr_df$y, Comparison = fdr_df$comparison, label = fdr_df$label)
            plt_ratio = plt_ratio + geom_text(data = ann_fdr_assoc, label = '**', vjust = -1.5)
            # annotation of mean change
            change_df = data.frame(gene = c(1,1,1), ratio_beta = c(15, 15, 13), Comparison = rep('Not different', 3), label = c('Known SNP', 'New SNP in latest GWAS', 'Known SNP'))
            labels_change = c(paste0('Median ratio known SNP = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta))) & tmp$label == 'Known SNP']), 2)),
                paste0('Median ratio new SNP = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta))) & tmp$label != 'Known SNP']), 2)),
                paste0('Median ratio overall = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta)))]), 2)))
            plt_ratio = plt_ratio + geom_text(data = change_df, label = labels_change, hjust = 0)
        # manually combine the plots together
        combined_plot = ggarrange(plotlist = list(plt, plt_grp, plt_ratio), nrow = 3, ncol = 1, labels = "AUTO")
        pdf('figure_1.pdf', height = 19, width = 16)
        combined_plot
        dev.off()
        return('** plot is done!')
    }

    # function to draw figure S1
    plot_figure_S1 <- function(singleAssoc, ratios, data_path, ad_snps){
        # create right order for the snps
            tmp = ratios[which(ratios$test == 'ad_vs_ctr')]
            tmp$ratio_beta = tmp$beta_aligned / tmp$gwas_beta; tmp = tmp[order(tmp$ratio_beta),]; 
            tmp$gene[which(duplicated(tmp$gene) == TRUE)] = paste0(tmp$gene[which(duplicated(tmp$gene) == TRUE)], ' (2)')
            tmp$gene = factor(tmp$gene, levels = tmp$gene)
            ad_snps_new_info = ad_snps[, c('locus', 'new_known')]
            tmp = merge(tmp, ad_snps_new_info, by = 'locus')
        # plot 1 is the grouped barplot with effect-sizes
            # prepare data for grouped barplot
            grp_plot = list()
            for (i in 1:nrow(tmp)){ grp_plot[[(length(grp_plot) + 1)]] = data.frame(Gene = rep(tmp$gene[i], 2), Effect = c(tmp$gwas_beta[i], tmp$beta_aligned[i]), Study = c('GWAS', 'AD cases vs. Centenarians'), low = c(tmp$gwas_beta[i] - (1.96*tmp$se_gwas[i]), tmp$beta_aligned[i] - (1.96*tmp$se[i])), up = c(tmp$gwas_beta[i] + (1.96*tmp$se_gwas[i]), tmp$beta_aligned[i] + (1.96*tmp$se[i]))) }
            grp_plot = rbindlist(grp_plot)
            # add labels for new of known snp
            new_snps = all_freqs_plot$gene[which(all_freqs_plot$new_known == 'new')]
            grp_plot$label = ifelse(grp_plot$Gene %in% new_snps, 'New SNP in latest GWAS', 'Known SNP')
            # plot
            plt_grp = ggplot(grp_plot, aes(fill=Study, y=Effect, x=Gene)) + geom_bar(position="dodge", stat="identity") + coord_cartesian(ylim = c(-1, 1), clip = "on") + xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + geom_errorbar(aes(ymin = low, ymax = up), width = 0.2, color = 'grey60', position=position_dodge(.9)) + ylim(-1, 1)
            plt_grp = plt_grp + facet_grid(cols = vars(label), drop = FALSE, scales = "free")
        # plot 3 is the effect-size ratio
            tmp$Comparison_symbol = ''; tmp$Comparison_symbol[which(tmp$two_tail_p <= 0.05)] = 'x'
            tmp$Comparison = 'Not different'; tmp$Comparison[which(tmp$two_tail_p <= 0.05)] = 'Different from published'
            tmp$adj_p = p.adjust(tmp$p, 'fdr')
            # add labels for new of known snp
            tmp$label = ifelse(tmp$gene %in% new_snps, 'New SNP in latest GWAS', 'Known SNP')
            # find significant associations to be annotated
            star_df = data.frame(x = tmp$gene[which(tmp$p <= 0.05)], y = tmp$up_ci_ratio[which(tmp$p <= 0.05)], comparison = tmp$Comparison[which(tmp$p <= 0.05)], label = tmp$label[which(tmp$p <= 0.05)])
            fdr_df = data.frame(x = tmp$gene[which(tmp$adj_p <= 0.05)], y = tmp$up_ci_ratio[which(tmp$adj_p <= 0.05)], comparison = tmp$Comparison[which(tmp$adj_p <= 0.05)], label = tmp$label[which(tmp$adj_p <= 0.05)])
            #tmp$ratio_beta[which(tmp$gene == 'TREML2')] = NA; tmp$low_ci_ratio[which(tmp$gene == 'TREML2')] = NA; tmp$up_ci_ratio[which(tmp$gene == 'TREML2')] = NA
            # the actual plot
            plt_ratio = ggplot(tmp, aes(x = gene, y = ratio_beta, fill = Comparison)) + geom_bar(stat = 'identity') + geom_errorbar(aes(ymin = low_ci_ratio, ymax = up_ci_ratio), width = 0.4, color = 'grey60') + geom_hline(yintercept = 1, linetype = 'dashed', col = 'red') +
                theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + xlab('') + ylab('Effect-size Ratio')
            plt_ratio = plt_ratio + facet_grid(cols = vars(label), drop = FALSE, scales = "free")
            # annotation of significant associations
            ann_sign_assoc = data.frame(gene = star_df$x, ratio_beta = star_df$y, Comparison = star_df$comparison, label = star_df$label)
            plt_ratio = plt_ratio + geom_text(data = ann_sign_assoc, label = '*')
            # annotation of fdr significant associations
            ann_fdr_assoc = data.frame(gene = fdr_df$x, ratio_beta = fdr_df$y, Comparison = fdr_df$comparison, label = fdr_df$label)
            plt_ratio = plt_ratio + geom_text(data = ann_fdr_assoc, label = '**', vjust = -1.5)
            # annotation of mean change
            change_df = data.frame(gene = c(1,1,1), ratio_beta = c(15, 15, 13), Comparison = rep('Not different', 3), label = c('Known SNP', 'New SNP in latest GWAS', 'Known SNP'))
            labels_change = c(paste0('Median ratio known SNP = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta))) & tmp$label == 'Known SNP']), 2)),
                paste0('Median ratio new SNP = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta))) & tmp$label != 'Known SNP']), 2)),
                paste0('Median ratio overall = ', round(median(tmp$ratio_beta[(!(is.na(tmp$ratio_beta)))]), 2)))
            plt_ratio = plt_ratio + geom_text(data = change_df, label = labels_change, hjust = 0)
        # manually combine the plots together
        combined_plot = ggarrange(plotlist = list(plt_grp, plt_ratio), nrow = 2, ncol = 1, labels = "AUTO")
        pdf('figure_s1.pdf', height = 12, width = 16)
        combined_plot
        dev.off()
        return('** plot is done!')
    }

    # function to draw figure 2 along with the forest plot of associations
    figure_2 <- function(prs_allSNPs, prs_allSNPs_noAPOE, assoc_all_combined, pheno){
        # densities
            # subset of individuals of interest
            tmp = prs_allSNPs[which(prs_allSNPs$PHENO %in% c('Centenarian', "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Probable_AD", "Possible_AD", "AD_path")),]
            tmp_noAPOE = prs_allSNPs_noAPOE[which(prs_allSNPs_noAPOE$PHENO %in% c('Centenarian', "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD", "Probable_AD", "Possible_AD", "AD_path")),]
            # take those that passed qc
            tmp = tmp[which(tmp$IID %in% pheno$ID_GWAS),]
            tmp_noAPOE = tmp_noAPOE[which(tmp_noAPOE$IID %in% pheno$ID_GWAS),]
            # manage children to keep
            children = fread('samples_children_controls.txt', h=F, stringsAsFactors=F)
            children_tokeep = prs_allSNPs[which(prs_allSNPs$IID %in% children$V2 & !(prs_allSNPs$IID %in% tmp$IID)),]
            children_tokeep_noAPOE = prs_allSNPs_noAPOE[which(prs_allSNPs_noAPOE$IID %in% children$V2 & !(prs_allSNPs_noAPOE$IID %in% tmp$IID)),]
            # add children to other samples
            #tmp = rbind(tmp, children_tokeep); tmp_noAPOE = rbind(tmp_noAPOE, children_tokeep_noAPOE)
            # add a label for centenarians, controls and AD
            tmp$Class = 'AD cases'; tmp$Class[which(tmp$PHENO == 'Centenarian')] = 'Centenarians'; tmp$Class[which(tmp$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 'Controls';
            tmp$Class[which(tmp$PHENO == 'family_100plus')] = 'Children'
            tmp_noAPOE$Class = 'AD cases'; tmp_noAPOE$Class[which(tmp_noAPOE$PHENO == 'Centenarian')] = 'Centenarians'; tmp_noAPOE$Class[which(tmp_noAPOE$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 'Controls'
            tmp_noAPOE$Class[which(tmp_noAPOE$PHENO == 'family_100plus')] = 'Children'
            # plot
            density1 = ggplot(tmp, aes(x = PRS, fill = Class)) + geom_density(alpha = 0.5) + ylab('Density') + ggtitle('PRS including APOE SNPs')
            density2 = ggplot(tmp_noAPOE, aes(x = PRS, fill = Class)) + geom_density(alpha = 0.5) + ylab('Density') + ggtitle('PRS excluding APOE SNPs')
            # combine plots
            combined = ggarrange(plotlist = list(density1, density2), nrow = 1, ncol = 2, common.legend = TRUE, labels = "AUTO")
        
        # forest plot
            # reformat data
            tmp_apoe = assoc_all_combined[which(assoc_all_combined$label == 'apoe_included'),]
            tmp_noapoe = assoc_all_combined[which(assoc_all_combined$label == 'apoe_excluded'),]
            # reorder data with title rows
            tmp_apoe_info = data.frame(comparison = 'All SNPs including APOE (N=86)', N_cases = "", N_ctr = "")
            tmp_noapoe_info = data.frame(comparison = 'All SNPs excluding APOE (N=84)', N_cases = "", N_ctr = "")
            combined_data = rbind.fill(tmp_apoe_info, tmp_apoe, tmp_noapoe_info, tmp_noapoe)
            # rename columns
            colnames(combined_data) = c('Subgroup', '# of cases', '# of controls', 'Beta', 'SE', 'Z', 'P', 'OR', 'lowCI', 'upCI', 'label', 'P (FDR)')
            # restrict to some columns only
            tmp_forplot = combined_data[, c('Subgroup', '# of cases', '# of controls', 'OR', 'lowCI', 'upCI', 'SE', 'P (FDR)')]
            # add width of the plot
            tmp_forplot$'    ' <- paste(rep("   ", nrow(tmp_forplot)), collapse = " ")
            # add the combined OR column
            tmp_forplot$"OR (95% CI)" = paste0(round(tmp_forplot$OR, 2), ' (', round(tmp_forplot$lowCI, 2), ' to ', round(tmp_forplot$upCI, 2), ')')
            tmp_forplot$"OR (95% CI)"[which(tmp_forplot$"OR (95% CI)" == 'NA (NA to NA)')] = ""
            # adjust labels
            tmp_forplot$Subgroup[which(tmp_forplot$Subgroup == 'ad_ctr')] = 'AD vs. Controls'
            tmp_forplot$Subgroup[which(tmp_forplot$Subgroup == 'chc_ad')] = 'AD vs. Centenarians'
            tmp_forplot$Subgroup[which(tmp_forplot$Subgroup == 'chc_ctr')] = 'Controls vs. Centenarians'
            tmp_forplot$Subgroup[which(tmp_forplot$Subgroup == 'children_vs_ctr')] = 'Controls vs. Children'
            # adjust pvalues
            tmp_forplot$"P (FDR)" = signif(tmp_forplot$"P (FDR)", digits = 3)
            tmp_forplot$"P (FDR)"[is.na(tmp_forplot$"P (FDR)")] = ""
            # plot
            tm = forest_theme(base_size = 15,
                        # center align column 2, 3, 5, 6
                        core=list(fg_params=list(hjust=c(0.5), x=c(0.5))),
                        # Confidence interval point shape, line type/color/width
                        ci_pch = 16,
                        ci_col = "#762a83",
                        ci_lty = 1,
                        ci_lwd = 2,
                        ci_Theight = 0.2, # Set an T end at the end of CI 
                        # Reference line width/type/color
                        refline_lwd = 1,
                        refline_lty = "dashed",
                        refline_col = "grey20",
                        # Vertical line width/type/color
                        vertline_lwd = 1,
                        vertline_lty = "dashed",
                        vertline_col = "grey20",
                        # Change summary color for filling and borders
                        summary_fill = "#4575b4",
                        summary_col = "#4575b4",
                        # Footnote font size/face/color
                        footnote_cex = 0.6,
                        footnote_fontface = "italic",
                        footnote_col = "blue")
            forestplot = forest(tmp_forplot[, c(1:3, 9:10, 8)], est = tmp_forplot$OR, lower = tmp_forplot$lowCI, upper = tmp_forplot$upCI, sizes = tmp_forplot$SE, ref_line = 1, ci_column = 4, arrow_lab = c('', 'Risk'), xlim = c(0, 6), ticks_at = c(0, 2, 4, 6), theme = tm)
            #forestplot_edit = edit_plot(forestplot, row = c(1, 6), gp = gpar(fontface = "bold"))
        # combine with densities and make plot
        combined2 = ggarrange(plotlist=list(combined, forestplot), nrow = 2, ncol = 1, labels = c('', 'C'), heights = c(0.50, 0.50))
        pdf('figure_2.pdf', height = 10, width = 16)
        combined2
        dev.off()
        return('plot is done!')
    }

    # extract imputation quality of the AD-snps
    extractQuality <- function(ad_snps, data_path){
        # change data path as info files are in another directory
        quality_path = str_replace_all(data_path, 'imputed_unscrambled_genotypes', 'imputed_genotypes')
        # define list to put results in
        quality_list = list()
        # loop across snps to find info
        for (i in 1:nrow(ad_snps)){
            cat(paste0('** processing SNP ', i, '\r'))
            cmd = paste0('zgrep chr', ad_snps$chrom[i], ':', ad_snps$pos[i], ': ', quality_path, '/chr', ad_snps$chrom[i], '.info.gz')
            info = system(cmd, intern = T); info = data.frame(str_split_fixed(info, '\t', 13))
            colnames(info) = c('snpid', 'ref', 'alt', 'alt_frq', 'maf', 'avg_call', 'rsq', 'genotyped', 'loorsq', 'empR', 'empRsq', 'dose0', 'dose1')
            quality_list[[(length(quality_list) + 1)]] = info
        }
        quality_list = rbindlist(quality_list)
        return(quality_list)
    }

    # function to run power analysis in parallel
    power_function = function(i, assoc, lit, n_cases, snplist, type, sampling, do_sampling){
        snp = snplist[i]
        df_converge = data.frame()
        df_all = data.frame()
        if (!(snp %in% c('1:109345810', '6:41181270', '6:41161514'))){
            # get association OR: ad vs. chc and ad vs. ctr
            or_chc = exp(assoc$beta_aligned[which(assoc$locus == snp & assoc$test == 'ad_vs_chc')])
            or_ctr = exp(assoc$beta_aligned[which(assoc$locus == snp & assoc$test == 'ad_vs_ctr')])
            # also do this in a sampling framework
            if (do_sampling == TRUE){
                or_chc_sample = exp(rnorm(mean = assoc$beta_aligned[which(assoc$locus == snp & assoc$test == 'ad_vs_chc')], sd = assoc$se[which(assoc$locus == snp & assoc$test == 'ad_vs_chc')], n=sampling))
                or_ctr_sample = exp(rnorm(mean = assoc$beta_aligned[which(assoc$locus == snp & assoc$test == 'ad_vs_ctr')], sd = assoc$se[which(assoc$locus == snp & assoc$test == 'ad_vs_ctr')], n=sampling))
            } else {
                or_chc_sample = or_chc
                or_ctr_sample = or_ctr
            }
            # get other parameters
            gene = unique(assoc$gene[which(assoc$locus == snp)])
            maf_gwas = lit$maf[which(lit$locus == snp)]
            alpha_chosen = 0.05; power = 0.80
            cat(paste0('** working on SNP --> ', i, ': ', gene, '\n'))
            # define the step for increasing number of controls and the maximum size of controls (2*ncases) -- now using 8000 AD cases and 200 samples as a step
            step_controls = 200
            max_controls = n_cases * 2
            # initialize sample sizes of controls, centenarians, and the power
            n_controls = 0
            desired_power = 0.80
            pw = 0
            # main loop for controls
            while (pw <= desired_power && n_controls <= (max_controls - step_controls)){
                # increase size of controls
                n_controls = n_controls + step_controls
                # calculate power
                if (type == 'centenarians'){
                    pw_ctr = custom_power(N = c(n_cases + n_controls), Case.Rate = n_cases/(n_cases + n_controls), k = NULL, MAF = maf_gwas, OR = or_chc_sample, Alpha = alpha_chosen, True.Model=c("Additive"), Test.Model=c("Additive"))
                } else {
                    pw_ctr = custom_power(N = c(n_cases + n_controls), Case.Rate = n_cases/(n_cases + n_controls), k = NULL, MAF = maf_gwas, OR = or_ctr_sample, Alpha = alpha_chosen, True.Model=c("Additive"), Test.Model=c("Additive"))
                }
                pw_ctr$Gene = snp
                # update the power
                if (is.na(pw_ctr$Power_at_Alpha_0.05)){
                    pw = 0
                } else {
                    pw = mean(na.omit(pw_ctr$Power_at_Alpha_0.05))
                }
                # save things
                df_all = rbind(df_all, pw_ctr)
            }
            # at the end, save the final result only
            df_converge = pw_ctr
            cat(paste0('**** Power for ', gene, ' is ', pw, ' with ', n_controls, '\n'))
            if (type == 'centenarians'){
                write.table(df_all, paste0('snp_', snp, '_power_iterations_centenarians.txt'), quote=F, row.names=F, sep = "\t")
            } else {
                write.table(df_all, paste0('snp_', snp, '_power_iterations.txt'), quote=F, row.names=F, sep = "\t")
            }
        }
        return(df_converge)
    }

    # custom function of package genpowr to make it work faster and without errors
    custom_power = function(N = NULL, Case.Rate = NULL, k = NULL, MAF = NULL, OR = NULL, Alpha = 0.05, True.Model = "All", Test.Model = "All"){
        # conditions to stop the calculation
        if (is.null(N) == T) { stop("N, the total sample size, must be specified.") }
        if (is.null(k) == T & is.null(Case.Rate) == T) { stop("k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, must be specified.") }
        if (is.null(k) == F & is.null(Case.Rate) == F) { stop("Specify one of k, the number of controls per case, or Case.Rate, the proportion of cases in the study sample, not both.") }
        if (is.null(MAF) == T) { stop("MAF (minor allele frequency) must be specified.") }
        if (is.null(OR) == T) { stop("OR (detectable odds ratio) must be specified.") }
        if (sum(Case.Rate >= 1) > 0 | sum(Case.Rate <= 0) > 0) { stop("R2 must be greater than 0 and less than 1.") }
        if (sum(MAF >= 1) > 0 | sum(MAF <= 0) > 0) { stop("MAF must be greater than 0 and less than 1.") }
        if (sum(N <= 0) > 0) { stop("N must be greater than 0.") }
        if (sum(k <= 0) > 0) { stop("k must be greater than 0.") }
        if (sum(OR <= 0) > 0) { stop("OR must be greater than 0.") }
        if (sum(Alpha >= 1) > 0 | sum(Alpha <= 0) > 0) { stop("Alpha must be greater than 0 and less than 1.") }
        if (sum(!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))) > 0) { stop(paste("Invalid Test.Model:", paste(Test.Model[!(Test.Model %in% c("Dominant", "Recessive", "Additive", "2df", "All"))], collapse = ", "))) }
        if (sum(!(True.Model %in% c("Dominant", "Recessive", "Additive", "Additive", "All"))) > 0) { stop(paste("Invalid True.Model:", paste(True.Model[!(True.Model %in% c("Dominant", "Recessive", "Additive", "All"))], collapse = ", "))) }
        if ("All" %in% Test.Model) { Test.Model <- c("Dominant", "Recessive", "Additive", "2df") }
        if ("All" %in% True.Model) { True.Model <- c("Dominant", "Recessive", "Additive") }
        if (is.null(Case.Rate) == T) { Case.Rate = 1/(1 + k) }
        # start of calculation
        sample.size.tab <- expand.grid(N, Case.Rate)
        colnames(sample.size.tab) <- c("N", "Case.Rate")
        sample.size.tab$N_cases <- sample.size.tab$N * sample.size.tab$Case.Rate
        sample.size.tab$N_controls <- sample.size.tab$N - sample.size.tab$N_cases
        iter <- nrow(sample.size.tab)
        final.pow.tab <- NULL
        for (zz in 1:iter) {
            N <- sample.size.tab[zz, "N"]
            Case.Rate <- sample.size.tab[zz, "Case.Rate"]
            N_cases <- sample.size.tab[zz, "N_cases"]
            N_controls <- sample.size.tab[zz, "N_controls"]
            o.save.tab <- NULL
            for (o in OR) {
                #cat(paste0('** ', o, '\n'))
                m.save.tab <- NULL
                for (m in MAF) {
                    save.tab <- NULL
                #    cat(paste0('*** P_AA - P_AB - P_BB \n'))
                    P_AA <- (1 - m)^2
                    P_AB <- 2 * m * (1 - m)
                    P_BB <- m^2
                    a <- (o - 1)
                    b <- (P_AB + o * P_BB + Case.Rate - Case.Rate * o)
                    c <- -P_AB * Case.Rate
                #    cat(paste0('**** Quad roots \n'))
                    soln <- quad_roots(a, b, c)[2]
                    upper.lim <- min(soln, P_AB)
                    fa.2 <- function(x) { o - x * (P_AA - Case.Rate + x + ((o * x * P_BB)/(P_AB - x + o * x)))/((Case.Rate - x - ((o * x * P_BB)/(P_AB - x + o * x))) * (P_AB - x)) }
                #    cat(paste0('***** Add roots \n'))
                    add2.root <- zero_finder_nleqslv_custom(fa.2, veclength = 1, x.start.vals = 0.5 * soln, upper.lim = upper.lim)
                    P_AB_case_a2 <- add2.root/P_AB
                    prob_AB_case_a2 <- P_AB_case_a2 * P_AB
                    prob_AB_control_a2 <- (1 - P_AB_case_a2) * P_AB
                    prob_AA_case_a2 <- P_AA * prob_AB_case_a2/(o * prob_AB_control_a2 + prob_AB_case_a2)
                    prob_AA_control_a2 <- P_AA - prob_AA_case_a2
                    prob_BB_case_a2 <- (prob_AA_case_a2 * P_BB * o^2)/(prob_AA_case_a2 * o^2 + prob_AA_control_a2)
                    prob_BB_control_a2 <- P_BB - prob_BB_case_a2
                #    cat(paste0('****** Condition 1 \n'))
                    if (length(add2.root) > 1) {
                        qqcr <- numeric(0)
                        for (aqq in 1:length(add2.root)) {
                            qqt <- rbind(c(prob_AA_case_a2[aqq], prob_AB_case_a2[aqq], prob_BB_case_a2[aqq]), c(prob_AA_control_a2[aqq], P_AB - prob_AB_case_a2[aqq], prob_BB_control_a2[aqq]))
                            qqcr <- c(qqcr, apply(qqt, 1, sum)[1])
                        }
                        qqz <- which(qqcr - Case.Rate < 0.001)
                        if (length(qqz) > 1) { qqz <- which.min(qqcr - Case.Rate < 0.001) }
                        if (length(qqz) > 1) { stop("trouble selecting correct zero in additive function") }
                        prob_AA_case_a2 <- prob_AA_case_a2[qqz]
                        prob_AB_case_a2 <- prob_AB_case_a2[qqz]
                        prob_BB_case_a2 <- prob_BB_case_a2[qqz]
                        prob_AA_control_a2 <- prob_AA_control_a2[qqz]
                        prob_AB_case_a2 <- prob_AB_case_a2[qqz]
                        prob_BB_control_a2 <- prob_BB_control_a2[qqz]
                    }
                #    cat(paste0('******* Saving \n'))
                    add.tab2 <- data.frame(model = rep("Additive", 2), table = rbind(c(prob_AA_case_a2, prob_AB_case_a2, prob_BB_case_a2), c(prob_AA_control_a2, P_AB - prob_AB_case_a2, prob_BB_control_a2)))
                    save.tab <- rbind(save.tab, add.tab2)
                    m.save.tab <- rbind(m.save.tab, data.frame(True.Model = save.tab[, 1], MAF = m, OR = o, Disease.Status = rep(c("case", "control"), nrow(save.tab)/2), Geno.AA = save.tab[, 2], Geno.AB = save.tab[, 3], Geno.BB = save.tab[, 4]))
                }
                o.save.tab <- rbind(o.save.tab, m.save.tab)
            }
            power.tab <- NULL
            for (mod in Test.Model) {
                temp <- NULL
                for (j in seq(1, nrow(o.save.tab), 2)) {
                    t <- o.save.tab[j:(j + 1), c("Geno.AA", "Geno.AB", 
                    "Geno.BB")]
                    ll.alt <- calc.like(logistic.mles(t, model = mod), 
                    t, model = mod)
                    ll.null <- null.ll(t)
                    stat <- 2 * (as.numeric(ll.alt - ll.null))
                    if (mod == "2df") {
                    pow = 1 - pchisq(qchisq(1 - Alpha, df = 2, 
                        ncp = 0), df = 2, ncp = N * stat)
                    }
                    else {
                    pow = pnorm(sqrt(N * stat) - qnorm(1 - Alpha/2)) + 
                        pnorm(-sqrt(N * stat) - qnorm(1 - Alpha/2)) * 
                        1
                    }
                    temp <- rbind(temp, pow)
                }
                power.tab <- rbind(power.tab, data.frame(Test.Model = mod, 
                    o.save.tab[seq(1, nrow(o.save.tab), 2), 1:3], 
                    N, N_cases, N_controls, Case.Rate, temp))
            }
            colnames(power.tab) <- c("Test.Model", "True.Model", "MAF", "OR", "N_total", "N_cases", "N_controls", "Case.Rate", paste("Power_at_Alpha_", Alpha, sep = ""))
            final.pow.tab <- rbind(final.pow.tab, power.tab)
        }
        return(final.pow.tab)
    }

    # custom function of package genpowr to fix error
    zero_finder_nleqslv_custom = function (afun, veclength, tol = 0.4, x.start.vals = NULL, upper.lim = Inf) {
        conv <- reps <- reps2 <- 0
        res0 <- Inf
        while (any(res0 > upper.lim) | any(res0 == Inf)) {
            reps2 <- reps2 + 1
            if (reps2 > 100) 
                stop("cannot find a solution under upper.lim")
            while (conv == 0) {
                reps <- reps + 1
                if (is.null(x.start.vals)) {
                    x.start <- c(runif(veclength) * tol)
                } else {
                    if (reps + reps2 >= 250) {
                    x.start <- c(runif(veclength) * tol)
                    }
                    else if (reps + reps2 >= 150) {
                    x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, 
                        sd = x, n = 1)))
                    }
                    else if (reps + reps2 >= 50) {
                    x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, 
                        sd = x/10, n = 1)))
                    }
                    else {
                    x.start <- sapply(x.start.vals, function(x) abs(rnorm(mean = x, sd = x/100, n = 1)))
                    }
                }
                res <- nleqslv(x.start, afun)
                if (res$message %in% c("Function criterion near zero", "x-values within tolerance 'xtol'") & all(res$x > 0) & all(res$x < 1)) { conv <- 1 }
            }
            res0 <- res$x
            if (any(res0 > upper.lim)) 
                conv <- 0
        }
        return(res0)
    }

    # function to check children -- keep 1 child per family
    checkchildren <- function(tmp_prs, pheno_final_raw){
        # separate children from controls
        children = tmp_prs[which(tmp_prs$PHENO == 'family_100plus'),]
        no_children = tmp_prs[which(tmp_prs$PHENO != 'family_100plus'),]
        # add 100plus family id
        children = merge(children, pheno_final_raw, by.x = 'IID', by.y = 'ID_GWAS')
        children$family_100plus = str_split_fixed(children$ID_100plus, '_', 3)[, 1]
        children = children[sample(nrow(children)),]
        children_unique = children[!duplicated(children$family_100plus),]
        # re-join with controls
        newdf = rbind.fill(children_unique, no_children[which(no_children$PHENO != 'family_100plus'),])
        # now we need to calculate pcs
        fam = fread('/project/holstegelab/Share/gwas_array/TOPMED/data/Pop_stratification_IBD/data_common_20k_random_hg19.fam', h=F, stringsAsFactors=F)
        fam_sb = fam[which(fam$V2 %in% newdf$IID),]
        write.table(fam_sb[, c('V1', 'V2')], 'samples_children_controls.txt', quote=F, row.names=F, col.names=F)
        system('plink --bfile /project/holstegelab/Share/gwas_array/TOPMED/data/Pop_stratification_IBD/data_common_20k_random_hg19 --keep samples_children_controls.txt --pca --out pca_children_controls')
        # read pca components
        pcs = fread('pca_children_controls.eigenvec', h=F, stringsAsFactors=F)
        colnames(pcs) = c('FID', 'IID', paste0('PC', seq(1, 20)))
        pcs_sb = pcs[, c('IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
        tmp_prs = merge(newdf, pcs_sb, by = 'IID')
        return(tmp_prs)
    }

# MAIN ANALYSES
# 1. read variants from AD paper -- both clinical AD and proxy data
    ad_snps <- fread("AD_snps.txt", h=T, stringsAsFactors=F)
    ad_snps_proxy <- fread("Bellenguez_2022.txt", h=T, stringsAsFactors=F)
    # rename columns
    colnames(ad_snps) = c('rsid', 'chrom', 'pos', 'gene', 'locus', 'minor/major', 'or (95% ci)', 'p', 'info')
    # add odds ratio and confidence interval
    ad_snps$or = str_split_fixed(ad_snps$"or (95% ci)", ' ', 2)[, 1]
    ad_snps$ci = str_split_fixed(ad_snps$"or (95% ci)", ' ', 2)[, 2]
    ad_snps$ci = str_replace_all(ad_snps$ci, '\\(', ''); ad_snps$ci = str_replace_all(ad_snps$ci, '\\)', '')
    # add variable for known/new
    ad_snps$new_known = ifelse(ad_snps$"old/new" == 'New', 'new', 'known')
    # add locus
    ad_snps_proxy$locus = paste0(ad_snps_proxy$chrom, ':', ad_snps_proxy$pos)
    ad_snps$locus = paste0('chr', ad_snps$chrom, ':', ad_snps$pos)

# 2. set the path to the imputed genotypes
    data_path = '/project/holstegelab/Share/gwas_array/TOPMED/data/imputed_unscrambled_genotypes'

# 3. set the path to the phenotype data
    pheno <- read.table("/project/holstegelab/Share/gwas_array/mapping_files/20211027_phenotypes_ADC_samplesKept.txt", h=T)
    load('/project/holstegelab/Share/gwas_array/mapping_files/phenotypes_20211027_All.Rdata')

# 4. extract snps from imputed data -- we use stage I SNP effects, as they were estimated using only clinically diagnosed samples
    dosages_pheno <- function_extractAndPrepare(ad_snps, data_path, pheno_final_raw)
    dosages_pheno_proxy <- function_extractAndPrepare(ad_snps_proxy, data_path, pheno_final_raw)

# 5. exclude 1 SNP (chr4:993555:GAGTT:G_GAGTT) as it should not be there
    dosages_pheno$'chr4:993555:GAGTT:G_GAGTT' = NULL; dosages_pheno_proxy$'chr4:993555:GAGTT:G_GAGTT' = NULL
    dosages_pheno_noAPOE = dosages_pheno_proxy
    dosages_pheno_noAPOE$"chr19:44908822:C:T_C" = NULL; dosages_pheno_noAPOE$"chr19:44908684:T:C_T" = NULL
    # at the end, we have: dosages with and without APOE of AD SNPs

# 6. then calculate prs -- this time will use the ad_snps_proxy estimates
    prs_allSNPs <- function_PRS(dosages_pheno_proxy, ad_snps_proxy)
    prs_allSNPs_noAPOE <- function_PRS(dosages_pheno_noAPOE, ad_snps_proxy)

# 7. run associations
    assoc_all = function_testAssoc(prs_allSNPs, pheno, pheno_final_raw, 'no_pca')
    assoc_all_noAPOE = function_testAssoc(prs_allSNPs_noAPOE, pheno, pheno_final_raw, 'no_pca')
    # then add odds ratios and ci
    assoc_all$OR = exp(assoc_all$Estimate); assoc_all$lowci = exp(assoc_all$Estimate - (1.96 * assoc_all$"Std. Error")); assoc_all$upci = exp(assoc_all$Estimate + (1.96 * assoc_all$"Std. Error"))
    assoc_all_noAPOE$OR = exp(assoc_all_noAPOE$Estimate); assoc_all_noAPOE$lowci = exp(assoc_all_noAPOE$Estimate - (1.96 * assoc_all_noAPOE$"Std. Error")); assoc_all_noAPOE$upci = exp(assoc_all_noAPOE$Estimate + (1.96 * assoc_all_noAPOE$"Std. Error"))
    # add labels and put all associations together
    assoc_all$label = 'apoe_included'; assoc_all_noAPOE$label = 'apoe_excluded'
    assoc_all_combined = rbind(assoc_all, assoc_all_noAPOE)
    # correct for multiple tests
    assoc_all_combined$p_adjust = p.adjust(assoc_all_combined$"Pr(>|z|", 'fdr')

# 8. single-variant associations (ctr vs chc -- ad vs chc -- ad vs ctr)
    singleAssoc = function_singleVar(dosages_pheno, pheno)
    assoc_chc_ctr = singleAssoc[which(singleAssoc$test == 'ctr_vs_chc'),]
    assoc_chc_ctr$p_adjust = p.adjust(assoc_chc_ctr$p, 'fdr')
    assoc_chc_ctr[which(assoc_chc_ctr$p_adjust < 0.05),]

# 9. effect size ratio
    ratios = effect_size_ratio(singleAssoc, ad_snps)
    # add chromosome and position
    ratios$chrom = stringr::str_split_fixed(ratios$locus, ':', 2)[, 1]
    ratios$pos = stringr::str_split_fixed(ratios$locus, ':', 2)[, 2]
    # save the effect-size ratios
    write.table(ratios, 'Effect_size_ratios_EADB_ADsnps.txt', quote=F, row.names=F, sep="\t")
    # extract snps information and merge with ratios
    snps_info = ad_snps_proxy[, c('rsid', 'new_known', 'locus', 'maf')]
    ratios_ad_chc = ratios[which(ratios$test == 'ad_vs_chc'),]
    ratios_ad_ctr = ratios[which(ratios$test == 'ad_vs_ctr'),]
    ratios_ad_chc = merge(ratios_ad_chc, snps_info, by = 'locus')
    ratios_ad_ctr = merge(ratios_ad_ctr, snps_info, by = 'locus')
    # save the effect ratios of chc and ctr separately
    write.table(ratios_ad_chc, 'Effect_size_ratios_AD_vs_CHC_clinicalAD.txt', quote=F, row.names=F, sep="\t", dec = ',')
    write.table(ratios_ad_ctr, 'Effect_size_ratios_AD_vs_CTR_clinicalAD.txt', quote=F, row.names=F, sep="\t", dec = ',')

# 10. check ad-chc association: nominally significant, adjust p-values
    # separate dataset of ctr vs. chc
    ctr_chc_ratios = ratios[which(ratios$test == 'ctr_vs_chc'),]
    # run analysis for ad_chc
    checkAssociations(ratios_ad_chc, ad_snps)
    # run analysis for ad_ctr
    checkAssociations(ratios_ad_ctr, ad_snps)
    # run analysis for ad_chc
    checkAssociations(ctr_chc_ratios, ad_snps)
    # adjust pvalues
    ratios_ad_chc$p_adjust = p.adjust(ratios_ad_chc$p, 'fdr')
    ratios_ad_ctr$p_adjust = p.adjust(ratios_ad_ctr$p, 'fdr')
    # compare ad-chc and ad-ctr effect-size change for new and known SNPs, and all
    wilcox.test(x = ratios_ad_chc$ratio_beta, y = ratios_ad_ctr$ratio_beta)
    wilcox.test(x = ratios_ad_chc$ratio_beta[which(ratios_ad_chc$new_known == 'new')], y = ratios_ad_ctr$ratio_beta[which(ratios_ad_ctr$new_known == 'new')])
    wilcox.test(x = ratios_ad_chc$ratio_beta[which(ratios_ad_chc$new_known == 'known')], y = ratios_ad_ctr$ratio_beta[which(ratios_ad_ctr$new_known == 'known')])

# 11. extract frequencies
    res_frequencies = extract_frequencies(ratios)
    all_freqs = res_frequencies[[1]]; all_freqs_plot = res_frequencies[[2]]

# 12. extract significant associations (single-variants) to do gene-set enrichment analysis
    sign_snp = ratios[which(ratios$test == 'ad_vs_chc' & ratios$p <= 0.05),]
    write.table(ad_snps$rsid[which(ad_snps$locus %in% sign_snp$locus)], 'significant_snp_ad_chc.txt', quote=F, row.names=F, col.names=F)

# 13. how many normal controls is a centenarian worth
    # identify snps in the same direction of effect as the reference gwas
    same_direction = ratios[which((ratios$beta_aligned * ratios$gwas_beta) >0),]
    same_direction_chc = same_direction[which(same_direction$test == 'ad_vs_chc'),]
    same_direction_ctr = same_direction[which(same_direction$test == 'ad_vs_ctr'),]
    snplist = same_direction_chc$locus[which(same_direction_chc$locus %in% same_direction_ctr$locus)]
    # read literature for MAF
    lit = ad_snps_proxy[, c('locus', 'rsid', 'maf', 'gene')]
    # run the function using normal controls and centenarians
    res_normal = mclapply(1:length(snplist), power_function, assoc = ratios, lit = lit, n_cases = 8000, snplist = snplist, type = 'normal', sampling = 100, do_sampling = TRUE, mc.cores = 8)
    res_centen = mclapply(1:length(snplist), power_function, assoc = ratios, lit = lit, n_cases = 8000, snplist = snplist, type = 'centenarians', sampling = 100, do_sampling = TRUE, mc.cores = 8)
    # same without sampling
    set.seed(549385)
    res_normal_nosampl = mclapply(1:length(snplist), power_function, assoc = ratios, lit = lit, n_cases = 8000, snplist = snplist, type = 'normal', sampling = 100, do_sampling = FALSE, mc.cores = 8)
    res_centen_nosampl = mclapply(1:length(snplist), power_function, assoc = ratios, lit = lit, n_cases = 8000, snplist = snplist, type = 'centenarians', sampling = 100, do_sampling = FALSE, mc.cores = 8)
    # combine results
    res_normal_combined = rbindlist(res_normal)
    res_centen_combined = rbindlist(res_centen)
    res_normal_combined_nosampl = rbindlist(res_normal_nosampl)
    res_centen_combined_nosampl = rbindlist(res_centen_nosampl)
    # change column names and then merge
    res_normal_combined_nosampl = res_normal_combined_nosampl[, c('Gene', 'MAF', 'OR', 'N_total', 'N_cases', 'N_controls', 'Case.Rate', 'Power_at_Alpha_0.05')]; colnames(res_normal_combined_nosampl) = c('locus', 'maf', 'or_ctr', 'n_total_ctr', 'n_cases_ctr', 'n_ctr', 'case_rate_ctr', 'power_ctr') 
    res_centen_combined_nosampl = res_centen_combined_nosampl[, c('Gene', 'OR', 'N_total', 'N_cases', 'N_controls', 'Case.Rate', 'Power_at_Alpha_0.05')]; colnames(res_centen_combined_nosampl) = c('locus', 'or_chc', 'n_total_chc', 'n_cases_chc', 'n_chc', 'case_rate_chc', 'power_chc') 
    res_all_combined_nosampl = merge(res_normal_combined_nosampl, res_centen_combined_nosampl, by = 'locus')
    # add label for converged/not converged
    res_all_combined_nosampl$converged_info = NA; res_all_combined_nosampl$converged_info[which(res_all_combined_nosampl$n_chc < 16000 & res_all_combined_nosampl$n_ctr < 16000)] = 'both'
    res_all_combined_nosampl$converged_info[which(res_all_combined_nosampl$n_ctr == 16000 & res_all_combined_nosampl$n_chc < 16000)] = 'ctr_not_converged'
    res_all_combined_nosampl$converged_info[which(res_all_combined_nosampl$n_ctr < 16000 & res_all_combined_nosampl$n_chc == 16000)] = 'chc_not_converged'
    res_all_combined_nosampl$converged_info[is.na(res_all_combined_nosampl$converged_info)] = 'both_not_converged'
    table(res_all_combined_nosampl$converged_info)
    res_all_combined_nosampl$ratio = res_all_combined_nosampl$n_ctr / res_all_combined_nosampl$n_chc
    # look at ratio -- mean and median
    mean(res_all_combined_nosampl$ratio[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    median(res_all_combined_nosampl$ratio[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    # look at controls number -- mean, sd and median
    mean(res_all_combined_nosampl$n_ctr[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    median(res_all_combined_nosampl$n_ctr[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    sd(res_all_combined_nosampl$n_ctr[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    # look at centenarian number -- mean, sd and median
    mean(res_all_combined_nosampl$n_chc[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    median(res_all_combined_nosampl$n_chc[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    sd(res_all_combined_nosampl$n_chc[which(res_all_combined_nosampl$converged_info != 'both_not_converged')])
    # select snps for functional analysis
    gsea_input = res_all_combined_nosampl$locus[which(res_all_combined_nosampl$ratio >3)]
    write.table(ad_snps_proxy$rsid[which(ad_snps_proxy$locus %in% gsea_input)], 'input_for_gsea.txt', quote=F, row.names=F, col.names=F)
    # add some annotation like gene and maf
    res_all_combined_nosampl = merge(res_all_combined_nosampl, snps_info[, c('locus', 'rsid')], by = 'locus')
    # add description
    res_all_combined_nosampl$Description = ifelse(res_all_combined_nosampl$ratio >1, 'Centenarians value more', 'Centenarians value less')
    res_all_combined_nosampl$Description[which(res_all_combined_nosampl$ratio == 1)] = 'Same value (did not converge)'
    # add gene names and fix order for plot
    gene_info = ad_snps[, c('locus', 'gene')]
    gene_info$gene[which(duplicated(gene_info$gene) == TRUE)] = paste0(gene_info$gene[which(duplicated(gene_info$gene) == TRUE)], ' (2)')
    res_all_combined_nosampl = merge(res_all_combined_nosampl, gene_info, by = 'locus')
    # order by ratio and fix order
    res_all_combined_nosampl = res_all_combined_nosampl[order(res_all_combined_nosampl$ratio),]
    res_all_combined_nosampl$gene = factor(res_all_combined_nosampl$gene, levels = res_all_combined_nosampl$gene)
    # for statistics, exclude the snps that did not converge
    converged_all = res_all_combined_nosampl[which(res_all_combined_nosampl$n_ctr != 16000 | res_all_combined_nosampl$n_chc != 16000),]
    # add new snps info
    new_snps = ad_snps$locus[which(ad_snps$new_known == 'new')]   
    converged_all$Type = ifelse(converged_all$locus %in% new_snps, 'New SNP in latest GWAS', 'Known SNP')
    # plot figure 3
    pipi = ggplot(data = converged_all, aes(x = gene, y = log(ratio), fill = Description)) + geom_bar(stat = 'identity') + 
        xlab("") + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + geom_hline(yintercept = 1, linetype = 'dashed', col = 'red')
    pipi = pipi + facet_grid(cols = vars(Type), scales = 'free')    
    change_df = data.frame(gene = c(1,1,1), ratio = c(40, 37, 40), Type = c('Known SNP', 'Known SNP', 'New SNP in latest GWAS'), Description = rep('Centenarians value more', 3))
    labels_change = c(paste0('Mean ratio known SNP = ', round(mean(converged_all$ratio[which(converged_all$Type == 'Known SNP')]), 2)),
                paste0('Mean ratio new SNP = ', round(mean(converged_all$ratio[which(converged_all$Type == 'New SNP in latest GWAS')]), 2)),
                paste0('Mean ratio overall = ', round(mean(converged_all$ratio), 2)))
    pipi = pipi + geom_text(data = change_df, label = c(labels_change), hjust = 0)
    pdf('figure_3.pdf', width = 16, height = 9)
    pipi
    dev.off()
    # plot figure S2 as the raw number of samples including the SNP that did not converge
    # need to restructure data
    dataplot = data.frame()
    for (i in 1:nrow(res_all_combined_nosampl)){
        dataplot = rbind(dataplot, data.frame(gene = res_all_combined_nosampl$gene[i], n = c(res_all_combined_nosampl$n_ctr[i], res_all_combined_nosampl$n_chc[i]), Type = c('Normal Controls', 'Centenarians'), desc = res_all_combined_nosampl$Description[i]))
    }
    res_all_combined_nosampl$Description[which(res_all_combined_nosampl$n_chc == 16000 & res_all_combined_nosampl$n_ctr < 16000)] = 'Centenarians did not converge'
    res_all_combined_nosampl$Description[which(res_all_combined_nosampl$n_chc < 16000 & res_all_combined_nosampl$n_ctr == 16000)] = 'Normal Controls did not converge'
    pupu = ggplot(res_all_combined_nosampl, aes(x = n_ctr, y = n_chc, color = Description)) + geom_jitter(size = 4, alpha = 0.6) + geom_abline(slope = 1, linetype = 'dashed', col = 'grey40') + xlab('Number of Normal Controls') +
        ylab('Number of Centenarians') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_colour_viridis_d()
    pdf('figure_s2.pdf', height = 7, width = 9)
    pupu
    dev.off()
    # maybe a different plot?
    # need again to re-structure data
    tmp_df = data.frame()
    tmp_df2 = data.frame()
    index_order = 1
    index_order2 = 1
    for (i in 1:nrow(converged_all)){
        tmp_df = rbind(tmp_df, data.frame(index = c(index_order, index_order+1), gene = rep(as.character(converged_all$gene[i]), 2), Description = rep(converged_all$Description[i], 2), Type = rep(converged_all$Type[i], 2), n = c(converged_all$n_chc[i], converged_all$n_ctr[i]), Label = c('Centenarians', 'Normal Controls')))
        index_order = index_order + 2
        tmp_df2 = rbind(tmp_df2, data.frame(index = index_order2, gene = converged_all$gene[i], n_chc = converged_all$n_chc[i], n_ctr = converged_all$n_ctr[i], Type = converged_all$Type[i], Description = converged_all$Description[i]))
        index_order2 = index_order2 + 1
    }
    # plot 1
    ggplot(tmp_df, aes(x = gene, y = n, color = Label)) + geom_point(size = 3) + xlab('') + ylab('Number of individuals') +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top')
    # plot 2
    plt = ggplot(tmp_df2, aes(x = gene)) + labs(x = "", y = "Number of individuals") + 
        geom_segment(aes(x = gene, y = n_chc, xend = gene, yend = n_ctr), size = 1) + 
        geom_point(aes(y = n_chc, color = "Centenarians"), size = 4, shape = 15) + geom_point(aes(y = n_ctr, color = "Normal Controls"), size = 4, shape = 15) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'top') + scale_color_viridis(discrete = T, option = "D", end = 0.75)
    plt = plt + facet_grid(cols = vars(Type), scales = 'free')    
    pdf('figure_3_alternative.pdf', height = 9, width = 16)
    plt
    dev.off()
    # finally output the table
    write.table(res_all_combined_nosampl, 'table_s6.txt', quote = F, row.names = F, sep = "\t", dec = ',')

# PLOTS
# 1. figure 1 should be the maf, effect-size and ratio plots
    plot_figure1(singleAssoc, ratios, ad_snps, all_freqs_plot)

# 2. figure 2 should be the density plot
    figure_2(prs_allSNPs, prs_allSNPs_noAPOE, assoc_all_combined, pheno)

# TABLES
# 1. extract demographics
    # AD
    ad_samples = prs_allSNPs[which(prs_allSNPs$PHENO %in% c("Probable_AD", "Possible_AD", "AD_path")),]; ad_samples_info = pheno[which(pheno$ID_GWAS %in% ad_samples$IID),]
    mean(na.omit(ad_samples_info$age)); sd(na.omit(ad_samples_info$age)); table(ad_samples_info$sex)/nrow(ad_samples_info)
    # Controls
    ctr_samples = prs_allSNPs[which(prs_allSNPs$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]; ctr_samples_info = pheno[which(pheno$ID_GWAS %in% ctr_samples$IID),]
    mean(na.omit(ctr_samples_info$age)); sd(na.omit(ctr_samples_info$age)); table(ctr_samples_info$sex)/nrow(ctr_samples_info)
    # Centenarians
    chc_samples = prs_allSNPs[which(prs_allSNPs$PHENO %in% c("Centenarian")),]; chc_samples_info = pheno[which(pheno$ID_GWAS %in% chc_samples$IID),]
    mean(na.omit(chc_samples_info$age)); sd(na.omit(chc_samples_info$age)); table(chc_samples_info$sex)/nrow(chc_samples_info)
    # Children
    children = fread('pca_children_controls.eigenvec', h=F, stringsAsFactors=F)
    children_pheno = merge(children, pheno_final_raw, by.x = 'V2', by.y = 'ID_GWAS')
    children_only = children_pheno[which(children_pheno$diagnosis == 'family_100plus'),]
    mean(na.omit(children_only$age)); sd(na.omit(children_only$age)); table(children_only$sex)/nrow(children_only)

# 2. check imputation quality of the SNPs
    imputation_quality = extractQuality(ad_snps, data_path)
    imputation_quality$locus = str_replace_all(imputation_quality$locus, 'chr', '')
    imputation_quality_snps = merge(imputation_quality, ad_snps, by = 'locus')
    write.table(imputation_quality_snps, 'imputation_quality_values.txt', quote=F, row.names = F, sep = '\t')

# 3. make table S3 of all associations and ratios
    ratios_table_s3 = ratios[which(ratios$test == 'ad_vs_chc'),]
    ratios_table_s3$gene[which(duplicated(ratios_table_s3$gene) == TRUE)] = paste0(ratios_table_s3$gene[which(duplicated(ratios_table_s3$gene) == TRUE)], ' (2)')
    ratios_table_s3$p_fdr = p.adjust(ratios_table_s3$p, 'fdr')
    ad_snps_proxy$locus = paste(ad_snps_proxy$chrom, ad_snps_proxy$pos, sep = ":"); snps_info = ad_snps_proxy[, c('locus', 'rsid', 'maf')]
    ratios_annot = merge(ratios_table_s3, snps_info, by = 'locus')
    write.table(ratios_annot, 'table_s3.txt', quote=F, row.names=F, sep="\t", dec=',')

# WORKSPACE
    save.image('20220930_workspace_final.RData')
