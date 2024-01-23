#############################################################
# FINAL VERSION OF THE SCRIPT FOR THE NEW EFFECT-SIZE PAPER #
#############################################################

# LIBRARIES
    library(data.table)
    library(stringr)
    library(parallel)
    library(ggplot2)
    library(viridis)
    library(ggpubr)
    library(plyr)
    library(forestploter)
    library(grid)
    library(genpwr)
    library(nleqslv)
    library(ggpubr)
    library(RColorBrewer)
    library(showtext)

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
        # match real children only
        children_ids = pheno_final_raw$ID_GWAS[grep('_F', pheno_final_raw$ID_100plus)]
        children = children[which(children$IID %in% children_ids),]
        # keep only samples that passed qc and add pcs
        pcs = pheno[, c('ID_GWAS', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'sex')]
        prs = merge(prs, pcs, by.x = 'IID', by.y = 'ID_GWAS')
        # chc vs. controls
        tmp_prs = prs[which(prs$PHENO %in% c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO == 'Centenarian')] = 0; print(table(tmp_prs$STATUS))
        model_chc_ctr = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_chc_ctr = as.data.frame(t(summary(model_chc_ctr)$coefficients[2, ]))
        res_chc_ctr$comparison = 'chc_ctr'
        res_chc_ctr$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_chc_ctr$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # gender strat
        # males
        tmp_prs_males = tmp_prs[which(tmp_prs$sex == 'M'),]; tmp_prs_females = tmp_prs[which(tmp_prs$sex == 'F'),]
        model_chc_ctr_males = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_males, family = 'binomial')
        res_chc_ctr_males = as.data.frame(t(summary(model_chc_ctr_males)$coefficients[2, ]))
        res_chc_ctr_males$comparison = 'chc_ctr_males'
        res_chc_ctr_males$N_cases = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 1),]); res_chc_ctr_males$N_ctr = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 0),])
        # females
        model_chc_ctr_females = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_females, family = 'binomial')
        res_chc_ctr_females = as.data.frame(t(summary(model_chc_ctr_females)$coefficients[2, ]))
        res_chc_ctr_females$comparison = 'chc_ctr_females'
        res_chc_ctr_females$N_cases = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 1),]); res_chc_ctr_females$N_ctr = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 0),])
        
        # chc vs. ad
        tmp_prs = prs[which(prs$PHENO %in% c("Centenarian", "Probable_AD", "Possible_AD", "AD_path")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO == 'Centenarian')] = 0; print(table(tmp_prs$STATUS))
        model_chc_ad = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_chc_ad = as.data.frame(t(summary(model_chc_ad)$coefficients[2, ]))
        res_chc_ad$comparison = 'chc_ad'
        res_chc_ad$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_chc_ad$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # gender strat
        # males
        tmp_prs_males = tmp_prs[which(tmp_prs$sex == 'M'),]; tmp_prs_females = tmp_prs[which(tmp_prs$sex == 'F'),]
        model_chc_ad_males = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_males, family = 'binomial')
        res_chc_ad_males = as.data.frame(t(summary(model_chc_ad_males)$coefficients[2, ]))
        res_chc_ad_males$comparison = 'chc_ad_males'
        res_chc_ad_males$N_cases = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 1),]); res_chc_ad_males$N_ctr = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 0),])
        # females
        model_chc_ad_females = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_females, family = 'binomial')
        res_chc_ad_females = as.data.frame(t(summary(model_chc_ad_females)$coefficients[2, ]))
        res_chc_ad_females$comparison = 'chc_ad_females'
        res_chc_ad_females$N_cases = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 1),]); res_chc_ad_females$N_ctr = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 0),])
        
        # ad vs. ctr
        tmp_prs = prs[which(prs$PHENO %in% c("Control_LASA", "Control_other_twin", "Control_100plus", "Control_path", "SCD", "Probable_AD", "Possible_AD", "AD_path")),]
        tmp_prs$STATUS = 1; tmp_prs$STATUS[which(tmp_prs$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 0; print(table(tmp_prs$STATUS))
        model_ad_ctr = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs, family = 'binomial')
        res_ad_ctr = as.data.frame(t(summary(model_ad_ctr)$coefficients[2, ]))
        res_ad_ctr$comparison = 'ad_ctr'
        res_ad_ctr$N_cases = nrow(tmp_prs[which(tmp_prs$STATUS == 1),]); res_ad_ctr$N_ctr = nrow(tmp_prs[which(tmp_prs$STATUS == 0),])
        # gender strat
        # males
        tmp_prs_males = tmp_prs[which(tmp_prs$sex == 'M'),]; tmp_prs_females = tmp_prs[which(tmp_prs$sex == 'F'),]
        model_ad_ctr_males = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_males, family = 'binomial')
        res_ad_ctr_males = as.data.frame(t(summary(model_ad_ctr_males)$coefficients[2, ]))
        res_ad_ctr_males$comparison = 'ad_ctr_males'
        res_ad_ctr_males$N_cases = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 1),]); res_ad_ctr_males$N_ctr = nrow(tmp_prs_males[which(tmp_prs_males$STATUS == 0),])
        # females
        model_ad_ctr_females = glm(STATUS ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp_prs_females, family = 'binomial')
        res_ad_ctr_females = as.data.frame(t(summary(model_ad_ctr_females)$coefficients[2, ]))
        res_ad_ctr_females$comparison = 'ad_ctr_females'
        res_ad_ctr_females$N_cases = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 1),]); res_ad_ctr_females$N_ctr = nrow(tmp_prs_females[which(tmp_prs_females$STATUS == 0),])

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
        } else if (do_pca == 'all_children'){
            # here we would use all the children. Problem is that there are no PCs for 97 of them, so i need to generate them again
            children_clean = pcs_children(tmp_prs, 'no_pca')
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
        gender_strat = list(res_chc_ad_females, res_chc_ad_males, res_ad_ctr_females, res_ad_ctr_males, res_chc_ctr_females, res_chc_ctr_males)
        all_res = list(results, gender_strat)
        return(all_res)
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
            gwas_allele = str_split_fixed(ad_snps$"ea/oa"[which(ad_snps$locus == snp_locus)], '/', 2)[, 1]
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
        return(ratios_tmp)
    }

    # function to extract frequencies for the maf plot
    extract_frequencies <- function(ratios){
        snps_interest = str_split_fixed(unique(ratios$snp), '_', 2)[, 1]
        tmp = str_split_fixed(snps_interest, ':', 4)
        df = data.frame(snpid = snps_interest, chrom = tmp[, 1], position = tmp[, 2])
        # also need to write samples to look at
        pheno <- read.table("path/to/QC/GWAS/samples_kept.txt", h=T)
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
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep children_samples.txt --extract ', chr, '_snps.txt --freq --out ', chr, '_freq_off'))
            # read frequencies back in
            chc_freq = fread(paste0(chr, '_freq_chc.afreq'), h=T)
            ctr_freq = fread(paste0(chr, '_freq_ctr.afreq'), h=T)
            ad_freq = fread(paste0(chr, '_freq_ad.afreq'), h=T)
            off_freq = fread(paste0(chr, '_freq_off.afreq'), h=T)
            # combine frequencies
            df_freq = data.frame(snpid = chc_freq$ID, ref = chc_freq$REF, alt = chc_freq$ALT, n_allele_chc = chc_freq$OBS_CT, n_allele_ctr = ctr_freq$OBS_CT, n_allele_ad = ad_freq$OBS_CT, n_allele_off = off_freq$OBS_CT, freq_chc = chc_freq$ALT_FREQS, freq_ctr = ctr_freq$ALT_FREQS, freq_ad = ad_freq$ALT_FREQS, freq_off = off_freq$ALT_FREQS)
            df_freq_plot = rbind(chc_freq, ctr_freq, ad_freq, off_freq); df_freq_plot$Phenotype = c(rep('CHC', nrow(chc_freq)), rep('CTR', nrow(ctr_freq)), rep('AD', nrow(ad_freq)), rep('OFF', nrow(off_freq)))
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

    # function to extract frequencies for the maf plot
    extract_counts <- function(ratios){
        snps_interest = str_split_fixed(unique(ratios$snp), '_', 2)[, 1]
        tmp = str_split_fixed(snps_interest, ':', 4)
        df = data.frame(snpid = snps_interest, chrom = tmp[, 1], position = tmp[, 2])
        # main loop to calculate allele counts
        all_counts = list()
        all_counts_plot = list()
        for (chr in unique(df$chrom)){
            tmp_snps = df[which(df$chrom == chr),]
            write.table(tmp_snps$snpid, paste0(chr, '_snps.txt'), quote=F, row.names=F, col.names=F)
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep chc_samples.txt --extract ', chr, '_snps.txt --freq counts --out ', chr, '_freq_chc'))
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep ad_samples.txt --extract ', chr, '_snps.txt --freq counts --out ', chr, '_freq_ad'))
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep ctr_samples.txt --extract ', chr, '_snps.txt --freq counts --out ', chr, '_freq_ctr'))
            cmd = system(paste0("plink2 --pfile ", data_path, '/', chr, '.dose.unscrambled --keep children_samples.txt --extract ', chr, '_snps.txt --freq counts --out ', chr, '_freq_off'))
            # read frequencies back in
            chc_freq = fread(paste0(chr, '_freq_chc.acount'), h=T)
            ctr_freq = fread(paste0(chr, '_freq_ctr.acount'), h=T)
            ad_freq = fread(paste0(chr, '_freq_ad.acount'), h=T)
            off_freq = fread(paste0(chr, '_freq_off.acount'), h=T)
            # combine frequencies
            df_freq = data.frame(snpid = chc_freq$ID, ref = chc_freq$REF, alt = chc_freq$ALT, n_allele_chc = chc_freq$OBS_CT, n_allele_ctr = ctr_freq$OBS_CT, n_allele_ad = ad_freq$OBS_CT, n_allele_off = off_freq$OBS_CT, alt_count_chc = chc_freq$ALT_CTS, alt_count_ctr = ctr_freq$ALT_CTS, alt_count_ad = ad_freq$ALT_CTS, alt_count_off = off_freq$ALT_CTS)
            df_freq_plot = rbind(chc_freq, ctr_freq, ad_freq, off_freq); df_freq_plot$Phenotype = c(rep('CHC', nrow(chc_freq)), rep('CTR', nrow(ctr_freq)), rep('AD', nrow(ad_freq)), rep('OFF', nrow(off_freq)))
            # save
            all_counts[[(length(all_counts) + 1)]] = df_freq
            all_counts_plot[[(length(all_counts_plot) + 1)]] = df_freq_plot
        }
        all_counts = rbindlist(all_counts)
        tmp = str_split_fixed(all_counts$snpid, ':', 4)
        all_counts$locus = paste0(tmp[, 1], ':', tmp[, 2]); all_counts$locus = str_replace_all(all_counts$locus, 'chr', '')
        all_counts_plot = rbindlist(all_counts_plot)
        tmp2 = str_split_fixed(all_counts_plot$ID, ':', 4)
        all_counts_plot$locus = paste0(tmp2[, 1], ':', tmp2[, 2]); all_counts_plot$locus = str_replace_all(all_counts_plot$locus, 'chr', '')
        res = list(all_counts, all_counts_plot)
        return(res)
    }

    # function to extract hardcalls plot
    count_prot_dosages <- function(dosages_pheno, samples, ad_snps_info){
        # samples = 'chc_samples.txt'
        tmp_samples = fread(samples, h=F)
        tmp_dos = dosages_pheno[which(dosages_pheno$IID %in% tmp_samples$V1),]
        for (c in 1:ncol(tmp_dos)){
            if (length(grep('chr', colnames(tmp_dos)[c])) >0){
                tested_all = str_split_fixed(colnames(tmp_dos)[c], '_', 2)[, 2]
                locus_snp = paste(str_split_fixed(str_split_fixed(colnames(tmp_dos)[c], '_', 2)[, 1], ':', 4)[, 1:2], collapse=':')
                if (tested_all != ad_snps_info$prot_allele[which(ad_snps_info$locus == locus_snp)]){
                    tmp_dos[, c] = 2 - tmp_dos[, ..c]
                }
            }
        }
        # keep only snps
        snps_index = grep('chr', colnames(tmp_dos))
        tmp_dos = tmp_dos[, ..snps_index]
        # count the total number of variants in which there's at least 1 protective allele
        count <- apply(tmp_dos > 0.5, 1, sum)
        # count the total number of variants in which there're 2 protective alleles
        count_homo <- apply(tmp_dos > 1.5, 1, sum)
        # also add sum of number of protective alleles
        tmp_dos$sum_protective = rowSums(tmp_dos)
        tmp_dos$n_variants_with_protective = count
        tmp_dos$n_homo_protective = count_homo
        return(tmp_dos)
    }

    # function to draw figure 1 -- 3 plots: maf, effect-size and ratio
    figure1_custom = function(singleAssoc, ratios, ad_snps, all_freqs_plot, snps_info){
        # get positions for the plot first
        plot_pos = as.numeric(barplot(ratios$ratio_beta[which(ratios$test == 'ad_vs_chc')]))
        # order data the way we want: known snps first, old snps then, within each group order by beta ratio
        ratios_chc = ratios[which(ratios$test == 'ad_vs_chc'),]
        all_freqs_plot = merge(all_freqs_plot, snps_info, by = 'locus')
        all_freqs_plot = merge(all_freqs_plot, ratios_chc, by = 'locus')
        new_snps = all_freqs_plot[which(all_freqs_plot$new_known == 'new'),]
        known_snps = all_freqs_plot[which(all_freqs_plot$new_known == 'known'),]
        new_snps = new_snps[order(new_snps$ratio_beta),]
        known_snps = known_snps[order(known_snps$ratio_beta),]
        data_plot = rbind(known_snps, new_snps)
        # set up frequcny always to the minor allele
        data_plot$FREQ_TOPLOT = ifelse(data_plot$ALT_FREQS < 0.5, data_plot$ALT_FREQS, 1-data_plot$ALT_FREQS)
        data_plot$ALLELE_TOPLOT = ifelse(data_plot$ALT_FREQS < 0.5, data_plot$ALT, data_plot$REF)
        # assign risk or protection to minor alleles
        data_plot$snp_effect = NA
        for (i in 1:nrow(data_plot)){
            effect_allel_gwas = data_plot$gwas_allele[i]
            effect_gwas = as.numeric(data_plot$gwas_beta[i])
            if (data_plot$ALLELE_TOPLOT[i] == effect_allel_gwas){
            if (effect_gwas > 0){ data_plot$snp_effect[i] = 'Risk Allele' } else { data_plot$snp_effect[i] = 'Protective Allele' }
            } else {
            if (effect_gwas < 0){ data_plot$snp_effect[i] = 'Risk Allele' } else { data_plot$snp_effect[i] = 'Protective Allele' }
            }
        }
        # assign significance
        data_plot$Comparison = 'Not different'; data_plot$Comparison[which(data_plot$two_tail_p <= 0.05)] = 'Different from published'
        data_plot$adj_p = p.adjust(data_plot$p, 'fdr')
        # fix gene names
        data_plot$gene[which(data_plot$gene == 'SLC24A4')] = 'RIN3'
        data_plot$gene[duplicated(data_plot$gene)] = paste0(data_plot$gene[duplicated(data_plot$gene)], ' (2)')
        # fix label colors
        lab_color = data.frame(gene = as.character(), col = as.character(), locus = as.character())
        for (i in 1:nrow(data_plot)){
            if (!(data_plot$locus[i] %in% lab_color$locus)){ lab_color = rbind(lab_color, data.frame(gene = data_plot$gene[i], col = ifelse(data_plot$new_know[i] == 'new', 'navy', 'darkred'), locus = data_plot$locus[i])) }
        }
        # also get snps
        snplist = unique(data_plot$locus[order(data_plot$ratio_beta)])
        # set up grid for the plot
        pdf('Downloads/figure1_newCentenarians.pdf', height = 12, width = 18)
        sysfonts::font_add("Geneva", "/User/nicco/Desktop/Geneva.ttf")
        showtext::showtext_auto()
        par(mfrow=c(2, 1), family = 'Geneva')
        # plot 1 is the frequency plot
        par(mar = c(0, 6, 8, 3))
        plot(0, 0, xlim=c(0, ceiling(max(plot_pos))), xaxs = 'i', ylim = c(0, 0.5), bty = 'n', xaxt = 'none', xlab='', ylab='Minor Allele Frequency', pch = 16, col = 'white', cex.lab = 1.75)
        # grid
        for (x in plot_pos){ abline(v = x, lwd = 0.25, col = 'grey80') }
        for (y in seq(0, 0.5, length.out = 10)){ abline(h = y, lwd = 0.25, col = 'grey80') }
        # loop on snps
        counter = 1
        wid = 0.6
        color_palette = brewer.pal(n=3, "Pastel1")
        for (snp in snplist){
            # plot ctr
            points(x = plot_pos[counter], y = data_plot$FREQ_TOPLOT[which(data_plot$locus == snp & data_plot$Phenotype == 'CTR')], pch=15, col="deepskyblue3", cex=1.5, xpd=TRUE)
            # plot ad
            points(x = plot_pos[counter], y = data_plot$FREQ_TOPLOT[which(data_plot$locus == snp & data_plot$Phenotype == 'AD')], pch=16, col="red", cex=1.5, xpd=TRUE)
            # plot chc
            points(x = plot_pos[counter], y = data_plot$FREQ_TOPLOT[which(data_plot$locus == snp & data_plot$Phenotype == 'CHC')], pch=17, col="darkolivegreen4", cex=1.5, xpd=TRUE)
            # risk or protective square on top
            square_col = ifelse(data_plot$snp_effect[which(data_plot$locus == snp)][1] == 'Protective Allele', 'black', 'white')
            rect(xleft = plot_pos[counter] - wid, ybottom = 0.5, xright = plot_pos[counter] + wid, ytop = 0.525, col = square_col, xpd=T)
            counter = counter +1
        }
        legend(x = plot_pos[1], y = 0.7, legend = c('AD cases', 'Age-matched controls', 'Cognitively Healthy Centenarians'), pch = c(16, 15, 17), col = c('red', 'deepskyblue3', 'darkolivegreen4'), bty = 'n', xpd=T, ncol = 3, cex = 1.50)
        legend(x = plot_pos[1], y = 0.65, legend = c('Protective allele', 'Risk allele'), pch = c(15, 0), col = c('black', 'black'), bty = 'n', xpd=T, ncol = 2, cex = 1.5)
        legend(x = plot_pos[1], y = 0.6, legend = c('Different from published', 'Not different'), pch = c(15, 15), col = c(color_palette[1], color_palette[2]), bty = 'n', xpd=T, ncol = 2, cex = 1.5)
        segments(x0 = 0, y0 = -0.02, x1 = max(plot_pos), y1 = -0.02, lwd=1.5, col='grey10')
        text(x = -7, y = 0.675, labels = 'A', cex = 1, font = 2, xpd=T)
        # then plot ratio
        par(mar = c(8, 6, 0, 3), family = 'Geneva')
        plot(0, 0, xlim=c(0, ceiling(max(plot_pos))), xaxs = 'i', ylim = c(-10, 15), bty = 'n', xaxt = 'none', xlab='', ylab='Effect-size Ratio', pch = 16, col = 'white', cex.lab = 1.75)
        # grid
        for (x in plot_pos){ abline(v = x, lwd = 0.25, col = 'grey80') }
        for (y in seq(-15, 15, length.out = 10)){ abline(h = y, lwd = 0.25, col = 'grey80') }
        # loop on snps
        counter = 1
        wid = 0.6
        for (snp in snplist){
            if (snp != '6:41181270'){
            # plot ratio
            bar_col = ifelse(data_plot$Comparison[which(data_plot$locus == snp)] == 'Not different', color_palette[2], color_palette[1])[1]
            rect(xleft = plot_pos[counter]-wid, ybottom = 0, xright = plot_pos[counter]+wid, ytop = data_plot$ratio_beta[which(data_plot$locus == snp)][1], col = bar_col)
            segments(x0 = plot_pos[counter], y0 = data_plot$low_ci_ratio[which(data_plot$locus == snp)][1], x1 = plot_pos[counter], y1 = data_plot$up_ci_ratio[which(data_plot$locus == snp)][1], col = 'grey10', lwd=1)
            # add significance in case
            if (data_plot$adj_p[which(data_plot$locus == snp)][1] <= 0.05){
                text(x = plot_pos[counter], y = data_plot$up_ci_ratio[which(data_plot$locus == snp)][1] + data_plot$up_ci_ratio[which(data_plot$locus == snp)][1]*0.10, labels = '**', cex = 1)
            } else if (data_plot$p[which(data_plot$locus == snp)][1] <= 0.05){
                text(x = plot_pos[counter], y = data_plot$up_ci_ratio[which(data_plot$locus == snp)][1] + data_plot$up_ci_ratio[which(data_plot$locus == snp)][1]*0.10, labels = '*', cex = 1)
            }
            } else {
            text(x = plot_pos[counter], y = 0, labels = 'NA', cex = 0.60, xpd=TRUE)
            }
            # add gene name
            text(x = plot_pos[counter], y = -11, labels = data_plot$gene[which(data_plot$locus == snp)][1], adj = 1, srt = 45, xpd=T, cex = 0.90, col = as.character(lab_color$col[which(lab_color$locus == snp)]))
            # risk or protective square on top
            counter = counter +1
        }
        abline(h=1, col='red', lwd=2, lty=2)
        text(x = plot_pos[1], y = 14, labels = paste0('Median ratio: ', round(median(ratios_chc$ratio_beta), 2)), adj = 0)
        text(x = plot_pos[length(plot_pos)], y = -17, labels = 'Genes known', xpd=T, adj = 1, cex = 0.90, col = 'darkred')
        text(x = plot_pos[length(plot_pos)], y = -18, labels = 'Genes found in latest GWAS', xpd=T, adj = 1, cex = 0.90, col = 'navy')
        text(x = plot_pos[length(plot_pos)], y = -19, labels = '**: FDR<10%         *: P<0.05', xpd=T, adj = 1, cex = 0.90)
        text(x = -7, y = 15, labels = 'B', cex = 1, font = 2, xpd=T)
        dev.off()
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
        #children_tokeep = prs_allSNPs[which(prs_allSNPs$PHENO == 'family_100plus'),]
        #children_tokeep_noAPOE = prs_allSNPs_noAPOE[which(prs_allSNPs_noAPOE$PHENO == 'family_100plus'),]
        # add children to other samples
        #tmp = rbind(tmp, children_tokeep); tmp_noAPOE = rbind(tmp_noAPOE, children_tokeep_noAPOE)
        # add a label for centenarians, controls and AD
        tmp$Class = 'AD cases'; tmp$Class[which(tmp$PHENO == 'Centenarian')] = 'Centenarians'; tmp$Class[which(tmp$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 'Controls';
        #tmp$Class[which(tmp$PHENO == 'family_100plus')] = 'Offspring'
        tmp_noAPOE$Class = 'AD cases'; tmp_noAPOE$Class[which(tmp_noAPOE$PHENO == 'Centenarian')] = 'Centenarians'; tmp_noAPOE$Class[which(tmp_noAPOE$PHENO %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD"))] = 'Controls'
        #tmp_noAPOE$Class[which(tmp_noAPOE$PHENO == 'family_100plus')] = 'Offspring'
        tmp$Linetype = factor(ifelse(tmp$Class != 'Offspring', "1", "2"))
        tmp_noAPOE$Linetype = factor(ifelse(tmp_noAPOE$Class != 'Offspring', "1", "2"))
        # modify class name
        tmp$Class[which(tmp$Class == 'Controls')] = 'Age-matched controls'
        tmp$Class[which(tmp$Class == 'Centenarians')] = 'Cognitively Healthy Centenarians'
        tmp_noAPOE$Class[which(tmp_noAPOE$Class == 'Controls')] = 'Age-matched controls'
        tmp_noAPOE$Class[which(tmp_noAPOE$Class == 'Centenarians')] = 'Cognitively Healthy Centenarians'
        # plot
        density1 = ggplot(tmp, aes(x = PRS, fill = Class)) + geom_density(aes(linetype=Linetype, color = Class), size = 1) + ylab('Density') + ggtitle('PRS including APOE SNPs') + scale_fill_manual(values=c(alpha("red", 0.5), alpha("darkolivegreen3", 0.5), alpha('deepskyblue2', 0.5))) + scale_linetype_manual(values=c("solid"), guide = "none") + scale_color_manual(values = c("red", "darkolivegreen3", "deepskyblue2")) +  theme_minimal() + theme(text = element_text(family = "Geneva"))
        density2 = ggplot(tmp_noAPOE, aes(x = PRS, fill = Class)) + geom_density(aes(linetype=Linetype, color = Class), size = 1) + ylab('Density') + ggtitle('PRS excluding APOE SNPs') + scale_fill_manual(values=c(alpha("red", 0.5), alpha("darkolivegreen3", 0.5), alpha('deepskyblue2', 0.5))) + scale_linetype_manual(values=c("solid"), guide = "none") + scale_color_manual(values = c("red", "darkolivegreen3", "deepskyblue2", "black")) + theme_minimal() + theme(text = element_text(family = "Geneva"))
        # combine plots
        combined = ggarrange(plotlist = list(density1, density2), nrow = 1, ncol = 2, common.legend = TRUE, labels = "AUTO")
  
        # forest plot
        # reformat data
        tmp_apoe = assoc_all_combined[which(assoc_all_combined$label == 'apoe_included' & assoc_all_combined$comparison != 'children_vs_ctr'),]
        tmp_noapoe = assoc_all_combined[which(assoc_all_combined$label == 'apoe_excluded' & assoc_all_combined$comparison != 'children_vs_ctr'),]
        # reorder data with title rows
        tmp_apoe_info = data.frame(comparison = 'All SNPs including APOE (N=86)', N_cases = "", N_ctr = "")
        tmp_noapoe_info = data.frame(comparison = 'All SNPs excluding APOE (N=84)', N_cases = "", N_ctr = "")
        combined_data = rbind.fill(tmp_apoe_info, tmp_apoe, tmp_noapoe_info, tmp_noapoe)
        # rename columns
        colnames(combined_data) = c('Comparison', '# of cases', '# of controls', 'Beta', 'SE', 'Z', 'P', 'OR', 'lowCI', 'upCI', 'label', 'P (FDR)')
        # restrict to some columns only
        tmp_forplot = combined_data[, c('Comparison', '# of cases', '# of controls', 'OR', 'lowCI', 'upCI', 'SE', 'P (FDR)')]
        # add width of the plot
        tmp_forplot$'    ' <- paste(rep("   ", nrow(tmp_forplot)), collapse = " ")
        # add the combined OR column
        tmp_forplot$"OR (95% CI)" = paste0(round(tmp_forplot$OR, 2), ' (', round(tmp_forplot$lowCI, 2), ' to ', round(tmp_forplot$upCI, 2), ')')
        tmp_forplot$"OR (95% CI)"[which(tmp_forplot$"OR (95% CI)" == 'NA (NA to NA)')] = ""
        # adjust labels
        tmp_forplot$Comparison[which(tmp_forplot$Comparison == 'ad_ctr')] = 'AD cases vs. Age-matched controls'
        tmp_forplot$Comparison[which(tmp_forplot$Comparison == 'chc_ad')] = 'AD cases vs. Cognitively Healthy Centenarians      '
        tmp_forplot$Comparison[which(tmp_forplot$Comparison == 'chc_ctr')] = 'Age-matched controls vs. Cognitively Healthy Centenarians        '
        #tmp_forplot$Comparison[which(tmp_forplot$Comparison == 'children_vs_ctr')] = 'Controls vs. Offspring'
        # adjust pvalues
        tmp_forplot$"P (FDR)" = signif(as.numeric(tmp_forplot$"P (FDR)"), digits = 3)
        tmp_forplot$"P (FDR)"[is.na(tmp_forplot$"P (FDR)")] = ""
        # add indentations
        tmp_forplot$Comparison = ifelse(tmp_forplot$Comparison %in% c('All SNPs including APOE (N=86)', 'All SNPs excluding APOE (N=84)'), tmp_forplot$Comparison, paste0('    ', tmp_forplot$Comparison))
        tmp_forplot$`OR (95% CI)` = paste0('    ', tmp_forplot$`OR (95% CI)`)
        # plot
        tm = forest_theme(base_size = 13, base_family = "Geneva",
                    # text justification
                    core=list(fg_params=list(hjust = c(0), x=c(0))),
                    colhead=list(fg_params=list(hjust=c(0, -2, -2, 0.5, 0.55, 0.2), x=c(0, -2, -2, 0.5, 0.55, 0.2))),
                    # Confidence interval point shape, line type/color/width
                    ci_pch = 16,
                    ci_col = c('navy'),
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
                    summary_col = "#4575b4")
        forestplot = forest(tmp_forplot[, c(1:3, 9:10, 8)], est = tmp_forplot$OR, lower = tmp_forplot$lowCI, upper = tmp_forplot$upCI, sizes = 1, ref_line = 1, ci_column = 4, xlab='Odds Ratio on AD', arrow_length = 0, xlim = c(0, 6), ticks_at = c(0, 2, 4, 6), theme = tm)
        forestplot
        # combine with densities and make plot
        combined2 = ggarrange(plotlist=list(combined, forestplot), nrow = 2, ncol = 1, labels = c('', 'C'), heights = c(0.50, 0.50))
        pdf('Downloads/figure2_newCentenarians_last.pdf', height = 10, width = 16)
        font_add("Geneva", "/User/nicco/Desktop/Geneva.ttf")
        showtext_auto()
        combined2
        dev.off()
        return('plot is done!')
    }

    # function to draw figure 3
    figure_3 = function(res_all_combined_nosampl, ad_snps){
        color_palette = brewer.pal(4, 'Paired')
        # add new/known
        info = ad_snps[, c('locus', 'new_known')]
        info$locus = str_replace_all(info$locus, 'chr', '')
        res_all_combined_nosampl = merge(res_all_combined_nosampl, info, by='locus')
        # change names
        res_all_combined_nosampl$Description[which(res_all_combined_nosampl$Description == 'Centenarians value less than Age-matched controls')] = 'Cognitively Healthy Centenarians value less than Age-matched controls'
        res_all_combined_nosampl$Description[which(res_all_combined_nosampl$Description == 'Centenarians value more than Age-matched controls')] = 'Cognitively Healthy Centenarians value more than Age-matched controls'
        res_all_combined_nosampl$Description[which(res_all_combined_nosampl$Description == 'Same value (did not converge)')] = 'Did not converge'
        res_all_combined_nosampl$Description[which(res_all_combined_nosampl$Description == 'Did not converge' & res_all_combined_nosampl$converged_info == 'both')] = 'Same value'
        # sort
        res_all_combined_nosampl = res_all_combined_nosampl[order(res_all_combined_nosampl$ratio),]
        # plot pos
        pos = as.numeric(barplot(res_all_combined_nosampl$ratio))
        # background plot
        pdf('Downloads/figure_3_centenarians.pdf', height = 8, width = 16)
        # first figure is the actual number of cases and controls
        par(mfrow=c(1,2), mar=c(6, 6, 5, 0), family = 'Geneva')
        plot(res_all_combined_nosampl$n_chc, res_all_combined_nosampl$n_ctr, bty = 'n', ylab='Age-matched controls to reach 80% power', xlab='Cognitively Healthy Centenarians to reach 80% power', pch = 16, col = 'white', cex.lab = 1.50, xaxt='none', yaxt='none')
        # axes
        ticks = c(100, 500, 1000, 2500, 5000, 7500, 10000, 13000, 16000)
        axis(side=1, at=ticks, labels=rep("", length(ticks)), cex=1.25)
        for (x in ticks){ text(x = x, y = -700, adj = c(1, 1), labels = x, xpd=T, srt=45, cex=0.80) }
        axis(side=2, at=ticks, labels=rep("", length(ticks)), cex=1.25)
        for (x in ticks){ text(x = -700, y = x, adj = c(1, 0.5), labels = x, xpd=T, cex= 0.80) }
        # grid
        for (x in seq(0, 16000, length.out=10)){ abline(v = x, lwd=0.25, col='grey80') }
        for (y in seq(0, 16000, length.out=10)){ abline(h = y, lwd=0.25, col='grey80') }
        # diagnoal line
        abline(a=0, b=1, lty = 2, col='red', lwd=1.5)
        # points
        chc_more = res_all_combined_nosampl[which(res_all_combined_nosampl$Description == 'Cognitively Healthy Centenarians value more than Age-matched controls'),]
        points(chc_more$n_chc, chc_more$n_ctr, pch = 16, col = alpha(color_palette[1], 0.75), cex=2.25)
        chc_less = res_all_combined_nosampl[which(res_all_combined_nosampl$Description == 'Cognitively Healthy Centenarians value less than Age-matched controls'),]
        points(chc_less$n_chc, chc_less$n_ctr, pch = 16, col = alpha(color_palette[2], 0.75), cex=2.25)
        same = res_all_combined_nosampl[which(res_all_combined_nosampl$Description == 'Same value'),]
        points(same$n_chc, same$n_ctr, pch = 16, col = alpha(color_palette[3], 0.75), cex=2.25)
        nonconv = res_all_combined_nosampl[which(res_all_combined_nosampl$Description == 'Did not converge'),]
        points(nonconv$n_chc, nonconv$n_ctr, pch = 16, col = alpha(color_palette[4], 0.75), cex=2.25)
        # figure number
        text(x = -3000, y = 18700, labels = 'A', cex = 1, font = 2, xpd=T)
        legend(x = 0, y = 18600, legend = c('Cognitively Healthy Centenarians value more', 'Cognitively Healthy Centenarians value less', 'Same value', 'Did not converge'), pch = 21, pt.bg = c(color_palette), cex = 0.8, col = 'black', bty = 'n', bg = "white", xpd=T, ncol=2, pt.cex=1.8)
        
        # second figure is the ratio
        par(mar=c(6, 6, 5, 0), family = 'Geneva')
        plot(0, 0, ylim=c(0, ceiling(max(pos))), yaxt='n', yaxs = 'i', xlim = c(0, 45), bty = 'n', xaxt = 'none', ylab='', xlab='Centenarians compared to Age-matched controls', pch = 16, col = 'white', cex.lab = 1.50)
        # yaxis
        axis(side=1, at=seq(0, 45, 5), labels=seq(0, 45, 5), cex=1.25)
        # grid
        for (x in pos){ segments(x0 = 0, y0 = x, x1 = 45, y1 = , lwd = 0.25, col = 'grey80') }
        for (y in seq(0, 45, length.out = 10)){ abline(v = y, lwd = 0.25, col = 'grey80') }
        counter = 1
        wid = 0.6
        snplist = unique(res_all_combined_nosampl$locus)
        # check SLC24A4 to RIN3
        res_all_combined_nosampl$gene = str_replace_all(res_all_combined_nosampl$gene, 'SLC24A4', 'RIN3')
        for (snp in snplist){
            # plot ratio
            if (res_all_combined_nosampl$ratio[which(res_all_combined_nosampl$locus == snp)] <1){
                col = color_palette[2]
            } else if (res_all_combined_nosampl$ratio[which(res_all_combined_nosampl$locus == snp)] >1){
                col = color_palette[1]
            } else {
                col = ifelse(res_all_combined_nosampl$n_chc[which(res_all_combined_nosampl$locus == snp)] == 16000, color_palette[4], color_palette[3])
            }
            rect(xleft = 0, ybottom = pos[counter]-wid, xright = res_all_combined_nosampl$ratio[which(res_all_combined_nosampl$locus == snp)][1], ytop = pos[counter]+wid, col = col, xpd=T)
            lab_col = ifelse(res_all_combined_nosampl$new_known[which(res_all_combined_nosampl$locus == snp)] == 'new', 'navy', 'darkred')
            text(x = -0.5, y = pos[counter], adj = c(1, 0.5), col = lab_col, xpd=T, labels = res_all_combined_nosampl$gene[which(res_all_combined_nosampl$locus == snp)], cex=0.50)
            counter = counter + 1
        }
        abline(v=1, lty=2, lwd=1.5, col = 'red')
        abline(v=2, lty=2, lwd=1.5, col = 'blue')
        legend(x = 0, y = ceiling(max(pos)) + ceiling(max(pos))*0.12, legend = c('Cognitively Healthy Centenarians value more', 'Cognitively Healthy Centenarians value less', 'Same value', 'Did not converge'), pch = 22, pt.bg = c(color_palette), cex = 0.8, col = 'black', bty = 'n', bg = "white", xpd=T, ncol=2, pt.cex=1.8)
        text(x = 45, y = ceiling(max(pos)) + ceiling(max(pos))*0.09, labels = 'Genes known', xpd=T, adj = 1, cex = 0.70, col = 'darkred')
        text(x = 45, y = ceiling(max(pos)) + ceiling(max(pos))*0.07, labels = 'Genes found in latest GWAS', xpd=T, adj = 1, cex = 0.70, col = 'navy')
        # figure number
        text(x = -6, y = ceiling(max(pos)) + ceiling(max(pos))*0.12, labels = 'B', cex = 1, font = 2, xpd=T)
        dev.off()
    }

    # function to draw figure 4c -- the whole Figure 4 is assembled separately
    figure_4c = function(){
        full_annot = read.table('Downloads/RESULTS_snpXplorer/snp_annotation.txt', h=T, stringsAsFactors = F, sep="\t")
        previous_gset = fread('/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Downloads_shared/RESULTS_12334/Enrichment_results.txt', h=T, stringsAsFactors = F)
        last_gset = fread('/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Downloads_shared/RESULTS_12334/Enrichment_results_withIntersection.txt', h=T, stringsAsFactors = F)
        
        dim(previous_gset)
        dim(last_gset)
        table(last_gset$term_id %in% previous_gset$term_id)
        table(previous_gset$term_id %in% last_gset$term_id)
        # check the ones that are different
        last_gset[which(!(last_gset$term_id %in% previous_gset$term_id)),]  # no-one is significant
        previous_gset[which(!(previous_gset$term_id %in% last_gset$term_id)),]  # no-one is significant
        
        # isolate significant
        sign_previous = previous_gset[which(previous_gset$avgP <= 0.05),]
        sign_last = last_gset[which(last_gset$avgP <= 0.05),]
        dim(sign_previous)
        dim(sign_last)
        table(sign_previous$term_id %in% sign_last$term_id)   # one of the previous is not there anymore
        table(sign_last$term_id %in% sign_previous$term_id)   # 3 of the new ones were not in previous
        
        # merge significant in previous with all in the last
        tmp_last = last_gset
        colnames(tmp_last) = paste0(colnames(tmp_last), '_LAST')
        sign_previous_with_last = merge(sign_previous, tmp_last, by.x = 'term_id', by.y = 'term_id_LAST')  
        dim(sign_previous)
        dim(sign_previous_with_last)
        
        # plot logp
        ggplot(sign_previous_with_last, aes(x = log10P, y = log10P_LAST)) + geom_point(stat = 'identity')   # perfect correlation basically
        
        # clean data -- keep previous genesets and p and intersection from last
        clean_sign_info = sign_previous_with_last[, c('term_name', 'term_id', 'avgP', 'log10P', 'intersection_LAST')]
        # also interesting to add the cluster number here
        clusters_info = fread('/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Downloads_shared/RESULTS_12334/geneSet_enrichment_results_and_clusters.txt', h=T, stringsAsFactors = F)
        cluster_info_tmp = clusters_info[, c('term_id', 'cluster_in_dendro')]
        clean_sign_info = merge(clean_sign_info, cluster_info_tmp, by = 'term_id')
        
        # nice to make a link between snps and pathways -- read snps and do this
        snps_annot = fread('/Users/nicco/Library/Mobile Documents/com~apple~CloudDocs/Downloads_shared/RESULTS_12334/snp_annotation.txt', h=T, stringsAsFactors = F)
        snps_pathways = data.frame()
        for (i in 1:nrow(snps_annot)){
            # get genes associated with snp
            tmp_genes = unlist(strsplit(snps_annot$geneList[i], ','))
            # loop on genes
            for (g in tmp_genes){
            if (g != ""){
                tmp_pathways = clean_sign_info[grep(g, clean_sign_info$intersection_LAST)]
                tmp_pathways$gene = g
                tmp_pathways$snpid = snps_annot$ID[i]
                tmp_pathways$source_annot = snps_annot$source_finalGenes[i]
                tmp_pathways$all_genes = snps_annot$geneList[i]
                tmp_pathways = tmp_pathways[, c('snpid', 'gene', 'all_genes', 'source_annot', 'term_id', 'term_name', 'cluster_in_dendro')]
                snps_pathways = rbind(snps_pathways, tmp_pathways)
            }
            }
        }
        snps_pathways$cluster_in_dendro = factor(snps_pathways$cluster_in_dendro)
        snps_pathways = snps_pathways[order(snps_pathways$cluster_in_dendro),]
        snps_pathways$term_name = factor(snps_pathways$term_name, levels = unique(snps_pathways$term_name))
        label_colors = data.frame(term_name = as.character(), color = as.character())
        for (i in 1:nrow(snps_pathways)){
            if (!(snps_pathways$term_name[i] %in% label_colors$term_name)){
            label_colors = rbind(label_colors, data.frame(term_name = snps_pathways$term_name[i], col = ifelse(snps_pathways$cluster_in_dendro[i] == "1", 'navy', 'darkred')))    
            }
        }
        snps_pathways$ID = paste0(snps_pathways$snpid, ' ~ ', snps_pathways$gene)
        
        # some snps are missing -- why is that
        missing_snps = snps_annot[which(!(snps_annot$ID %in% snps_pathways$snpid)),]
        missing_info = data.frame()
        for (i in 1:nrow(missing_snps)){
            tmp_genes = unlist(strsplit(missing_snps$geneList[i], ','))
            for (g in tmp_genes){
            gene_info = last_gset[grep(g, last_gset$intersection), c('term_name', 'avgP')]
            gene_info$gene = g
            missing_info = rbind(missing_info, gene_info)
            }
        }
        head(missing_info[which(missing_info$gene == 'EPDR1'),])
        head(missing_info[which(missing_info$gene == 'ANKH'),])
        
        # how can we plot this?
        # at gene level
        ggplot(snps_pathways, aes(x = term_name, y = ID, fill = cluster_in_dendro)) + geom_tile(colour = 'black') + xlab('') + ylab('') + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, color = label_colors$col), legend.position = 'top')
        # or at SNP level
        tmp_ad_snps = ad_snps[, c('rsid', 'gene')]
        snps_pathways = merge(snps_pathways, tmp_ad_snps, by.x = 'snpid', by.y='rsid')
        colnames(snps_pathways)[which(colnames(snps_pathways) == 'cluster_in_dendro')] = 'Cluster'
        snps_pathways$"Cluster" = ifelse(snps_pathways$Cluster == 1, 'Immune activation/regulation', 'Endo-lysosomal/Clearance')
        pdf('figure_4c.pdf', height = 8, width = 14)
        ggplot(snps_pathways, aes(x = term_name, y = gene.y, fill = Cluster)) + geom_tile(colour = 'black') + xlab('') + ylab('') + theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = label_colors$col, size = 10), legend.position = 'top', text = element_text(family = "Geneva")) + scale_fill_manual(values = c('darkred', 'navy'))
        dev.off()
        
        # cluster 1 -- immune
        # list all pathways of this cluster
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "1"), c("term_name")])
        # snps involved in immune response
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "1" & snps_pathways$term_name %in% c("activation of immune response", "immune system process", "regulation of immune system process", "positive regulation of immune system process", "immune response", "regulation of immune response")), c("gene", "snpid")])
        # snps involved in laukocyte things
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "1" & snps_pathways$term_name %in% c("leukocyte activation", "leukocyte differentiation")), c("gene", "snpid")])
        # snps involved in macrophage things
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "1" & snps_pathways$term_name %in% c("macrophage activation involved in immune response", "macrophage activation")), c("gene", "snpid")])
        # snps involved in neuroinflamamtory things
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "1" & snps_pathways$term_name %in% c("positive regulation of inflammatory response", "neuroinflammatory response")), c("gene", "snpid")])
        
        # cluster 2 -- endocytosis
        # list all pathways of this cluster
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "2"), c("term_name")])
        # snps involved in endocytosis and phagocytosis
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "2" & snps_pathways$term_name %in% c("endocytosis", "phagocytosis, engulfment", "membrane invagination", "plasma membrane invagination", "regulation of phagocytosis, engulfment", "positive regulation of phagocytosis, engulfment", "regulation of membrane invagination", "positive regulation of membrane invagination")), c("gene", "snpid")])
        # snps involved in interleukine-6 things
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "2" & snps_pathways$term_name %in% c("interleukin-6 production", "regulation of interleukin-6 production", "positive regulation of interleukin-6 production")), c("gene", "snpid")])
        # snps involved in clearance
        unique(snps_pathways[which(snps_pathways$cluster_in_dendro == "2" & snps_pathways$term_name %in% c("amyloid-beta clearance by cellular catabolic process", "learning or memory", "neuropeptide processing")), c("gene", "snpid")])
        
        # completely different -- check how many of the novel SNPs are in the expected direction and test it
        novel_snps_ad_chc = ratios[which(ratios$test == 'ad_vs_chc'),]
        novel_snps_ad_chc = merge(novel_snps_ad_chc, snps_info, by = 'locus')
        novel_snps_ad_chc_new = novel_snps_ad_chc[which(novel_snps_ad_chc$new_known == 'new'),]
        nrow(novel_snps_ad_chc_new[which(novel_snps_ad_chc_new$ratio_beta >0),])/nrow(novel_snps_ad_chc_new)*100
        binom.test(x = nrow(novel_snps_ad_chc_new[which(novel_snps_ad_chc_new$ratio_beta >0),]), n = nrow(novel_snps_ad_chc_new), p = 0.5)
        
        # check overall significance
        singleAssoc_ad_chc = singleAssoc[which(singleAssoc$test == 'ad_vs_chc'),]
        singleAssoc_ad_ctr = singleAssoc[which(singleAssoc$test == 'ad_vs_ctr'),]
        singleAssoc_ctr_chc = singleAssoc[which(singleAssoc$test == 'ctr_vs_chc'),]
        singleAssoc_ad_chc$p_fdr = p.adjust(singleAssoc_ad_chc$p, 'fdr')
        singleAssoc_ad_ctr$p_fdr = p.adjust(singleAssoc_ad_ctr$p, 'fdr')
        singleAssoc_ctr_chc$p_fdr = p.adjust(singleAssoc_ctr_chc$p, 'fdr')
        dim(singleAssoc_ad_chc[which(singleAssoc_ad_chc$p_fdr <= 0.05),])
        dim(singleAssoc_ad_ctr[which(singleAssoc_ad_ctr$p_fdr <= 0.05),])
        dim(singleAssoc_ctr_chc[which(singleAssoc_ctr_chc$p_fdr <= 0.05),])
        sign_ad_chc = singleAssoc_ad_chc[which(singleAssoc_ad_chc$p_fdr <= 0.05), "snp"]
        sign_ad_ctr = singleAssoc_ad_ctr[which(singleAssoc_ad_ctr$p_fdr <= 0.05), "snp"]
        sign_ctr_chc = singleAssoc_ctr_chc[which(singleAssoc_ctr_chc$p_fdr <= 0.05), "snp"]
        intersect(sign_ad_chc, sign_ad_ctr)
        as.character(setdiff(sign_ad_chc, sign_ad_ctr)$snp)
        as.character(setdiff(sign_ad_ctr, sign_ad_chc)$snp)
        # check a couple
        # hla
        singleAssoc[which(singleAssoc$snp == 'chr6:32615322:A:G_A'),]   # nominal significant in ad vs. chc -- HLA
        singleAssoc[which(singleAssoc$snp == 'chr1:109345810:T:C_T'),]  # no significant -- SORT1
        singleAssoc[which(singleAssoc$snp == 'chr6:41161514:C:T_C'),]   # no significant -- TREM2
        singleAssoc[which(singleAssoc$snp == 'chr11:86157598:T:C_T'),]  # picalm nominal significant in ad vs. chc -- PICALM
        singleAssoc[which(singleAssoc$snp == 'chr15:63277703:C:T_C'),]  # no significant -- APH1B
        singleAssoc[which(singleAssoc$snp == 'chr16:81739398:G:A_G'),]  # borderline -- PLCG2
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

# MAIN ANALYSES
# 1. read variants from AD paper -- both clinical AD and proxy data
    ad_snps <- fread("AD_snps_clinicalAD.txt", h=T, stringsAsFactors=F)
    ad_snps_proxy <- fread("AD_snps_proxyAD.txt", h=T, stringsAsFactors=F)
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
    data_path = 'path/to/imputed/genotypes'

# 3. set the path to the phenotype data
    pheno <- read.table("path/to/QC/GWAS/samples_kept.txt", h=T)
    load('path/to/phenotype/file')

# 4. extract snps from imputed data -- we use stage I SNP effects, as they were estimated using only clinically diagnosed samples
    dosages_pheno <- function_extractAndPrepare(ad_snps, data_path, pheno_final_raw)
    dosages_pheno_proxy <- function_extractAndPrepare(ad_snps_proxy, data_path, pheno_final_raw)

# 5. exclude 1 SNP (chr4:993555:GAGTT:G_GAGTT) as it is multiallelic -- exclude the allele that is not associated with AD
    dosages_pheno$'chr4:993555:GAGTT:G_GAGTT' = NULL; dosages_pheno_proxy$'chr4:993555:GAGTT:G_GAGTT' = NULL
    dosages_pheno_noAPOE = dosages_pheno_proxy
    dosages_pheno_noAPOE$"chr19:44908822:C:T_C" = NULL; dosages_pheno_noAPOE$"chr19:44908684:T:C_T" = NULL
    # at the end, we have: dosages with and without APOE of AD SNPs

# 6. then calculate prs -- this time will use the ad_snps_proxy estimates
    prs_allSNPs <- function_PRS(dosages_pheno_proxy, ad_snps_proxy)
    prs_allSNPs_noAPOE <- function_PRS(dosages_pheno_noAPOE, ad_snps_proxy)

# 7. run associations -- association with children is still included despite not bein in the manuscript anymore
    resAssoc = function_testAssoc(prs_allSNPs, pheno, pheno_final_raw, 'all_children')
    assoc_all = resAssoc[[1]]; assoc_gender = rbindlist(resAssoc[[2]])
    resAssoc_noAPOE = function_testAssoc(prs_allSNPs_noAPOE, pheno, pheno_final_raw, 'all_children')
    assoc_all_noAPOE = resAssoc_noAPOE[[1]]; assoc_gender_noAPOE = rbindlist(resAssoc_noAPOE[[2]])
    # odds ratio and confidence intervals for the gender
    assoc_gender$or = exp(assoc_gender$Estimate); assoc_gender$low_ci = exp(assoc_gender$Estimate - 1.96*assoc_gender$"Std. Error"); assoc_gender$up_ci = exp(assoc_gender$Estimate + 1.96*assoc_gender$"Std. Error")
    assoc_gender_noAPOE$or = exp(assoc_gender_noAPOE$Estimate); assoc_gender_noAPOE$low_ci = exp(assoc_gender_noAPOE$Estimate - 1.96*assoc_gender_noAPOE$"Std. Error"); assoc_gender_noAPOE$up_ci = exp(assoc_gender_noAPOE$Estimate + 1.96*assoc_gender_noAPOE$"Std. Error")
    # same but keeping unrelated children instead of all of them
    assoc_all_unrelated = function_testAssoc(prs_allSNPs, pheno, pheno_final_raw, 'no_pca'); assoc_all_unrelated = assoc_all_unrelated[which(assoc_all_unrelated$comparison == 'children_vs_ctr'),]
    assoc_all_unrelated_noAPOE = function_testAssoc(prs_allSNPs_noAPOE, pheno, pheno_final_raw, 'no_pca'); assoc_all_unrelated_noAPOE = assoc_all_unrelated_noAPOE[which(assoc_all_unrelated_noAPOE$comparison == 'children_vs_ctr'),]
    # combine
    assoc_all = rbind(assoc_all, assoc_all_unrelated)
    assoc_all_noAPOE = rbind(assoc_all_noAPOE, assoc_all_unrelated_noAPOE)
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
    corrected_res_ad_chc = checkAssociations(ratios_ad_chc, ad_snps)
    # run analysis for ad_ctr
    corrected_res_ad_ctr = checkAssociations(ratios_ad_ctr, ad_snps)
    # run analysis for ad_chc
    corrected_res_ctr_chc = checkAssociations(ctr_chc_ratios, ad_snps)

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
    labels_change = c(paste0('Mean ratio known SNP = ', round(mean(converged_all$ratio[which(converged_all$Type == 'Known SNP')]), 2)), paste0('Mean ratio new SNP = ', round(mean(converged_all$ratio), 2)), paste0('Mean ratio overall = ', round(mean(converged_all$ratio[which(converged_all$Type == 'New SNP in latest GWAS')]), 2)))
    pipi = pipi + geom_text(data = change_df, label = c(labels_change), hjust = 0)
    pdf('figure_3.pdf', width = 16, height = 9)
    pipi
    dev.off()
    # plot figure S2 as the raw number of samples including the SNP that did not converge
    # need to restructure data
    dataplot = data.frame()
    for (i in 1:nrow(res_all_combined_nosampl)){ dataplot = rbind(dataplot, data.frame(gene = res_all_combined_nosampl$gene[i], n = c(res_all_combined_nosampl$n_ctr[i], res_all_combined_nosampl$n_chc[i]), Type = c('Normal Controls', 'Centenarians'), desc = res_all_combined_nosampl$Description[i])) }
    res_all_combined_nosampl$Description[which(res_all_combined_nosampl$n_chc == 16000 & res_all_combined_nosampl$n_ctr < 16000)] = 'Centenarians did not converge'
    res_all_combined_nosampl$Description[which(res_all_combined_nosampl$n_chc < 16000 & res_all_combined_nosampl$n_ctr == 16000)] = 'Normal Controls did not converge'
    pupu = ggplot(res_all_combined_nosampl, aes(x = n_ctr, y = n_chc, color = Description)) + geom_jitter(size = 4, alpha = 0.6) + geom_abline(slope = 1, linetype = 'dashed', col = 'grey40') + xlab('Number of Normal Controls') + ylab('Number of Centenarians') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_colour_viridis_d()
    pdf('figure_s2.pdf', height = 7, width = 9)
    pupu
    dev.off()
    # finally output the table
    write.table(res_all_combined_nosampl, 'table_s6.txt', quote = F, row.names = F, sep = "\t", dec = ',')

# 14. last thing: estimate how many protective allele each centenarian has compared to controls
    # load environment
    load('20230505_workspace_final.RData')
    write.table(children$IID, 'children_samples.txt', quote=F, row.names=F, col.names=F)
    # extract allele counts
    res_counts = extract_counts(ratios)
    allele_counts = res_counts[[1]]; allele_counts_plot = res_counts[[2]]
    # identify protective alleles
    allele_counts$locus = paste0('chr', allele_counts$locus)
    allele_counts$prot_allele_chc = NA; allele_counts$prot_allele_ctr = NA; allele_counts$prot_allele_ad = NA; allele_counts$prot_allele_off = NA;
    ad_snps_info = ad_snps[, c('minor/major', 'locus', 'or')]
    ad_snps_info$a1 = str_split_fixed(ad_snps_info$"minor/major", '/', 2)[, 1]
    ad_snps_info$a2 = str_split_fixed(ad_snps_info$"minor/major", '/', 2)[, 2]
    ad_snps_info$prot_allele = ifelse(ad_snps_info$or > 1, ad_snps_info$a2, ad_snps_info$a1)
    allele_counts = merge(allele_counts, ad_snps_info, by = 'locus')
    for (i in 1:nrow(allele_counts)){
        if (allele_counts$a1[i] == allele_counts$prot_allele[i]){
            allele_counts$prot_allele_chc[i] = allele_counts$alt_count_chc[i]
            allele_counts$prot_allele_ctr[i] = allele_counts$alt_count_ctr[i]
            allele_counts$prot_allele_ad[i] = allele_counts$alt_count_ad[i]
            allele_counts$prot_allele_off[i] = allele_counts$alt_count_off[i]
        } else {
            allele_counts$prot_allele_chc[i] = allele_counts$n_allele_chc[i] - allele_counts$alt_count_chc[i]
            allele_counts$prot_allele_ctr[i] = allele_counts$n_allele_ctr[i] - allele_counts$alt_count_ctr[i]
            allele_counts$prot_allele_ad[i] = allele_counts$n_allele_ad[i] - allele_counts$alt_count_ad[i]
            allele_counts$prot_allele_off[i] = allele_counts$n_allele_off[i] - allele_counts$alt_count_off[i]
        }
    }
    # maybe we should do it at the dosage level: question is how many protective alleles has each sample
    chc_sum_protective = count_prot_dosages(dosages_pheno, 'chc_samples.txt', ad_snps_info)
    ctr_sum_protective = count_prot_dosages(dosages_pheno, 'ctr_samples.txt', ad_snps_info)
    ad_sum_protective = count_prot_dosages(dosages_pheno, 'ad_samples.txt', ad_snps_info)
    off_sum_protective = count_prot_dosages(dosages_pheno, 'children_samples.txt', ad_snps_info)
    summary(chc_sum_protective$sum_protective)
    summary(off_sum_protective$sum_protective)
    summary(ctr_sum_protective$sum_protective)
    summary(ad_sum_protective$sum_protective)
    summary(chc_sum_protective$n_variants_with_protective)
    summary(off_sum_protective$n_variants_with_protective)
    summary(ctr_sum_protective$n_variants_with_protective)
    summary(ad_sum_protective$n_variants_with_protective)
    summary(chc_sum_protective$n_homo_protective)
    summary(off_sum_protective$n_homo_protective)
    summary(ctr_sum_protective$n_homo_protective)
    summary(ad_sum_protective$n_homo_protective)
    # combine
    df = data.frame(Value = c(chc_sum_protective$sum_protective, off_sum_protective$sum_protective, ctr_sum_protective$sum_protective, ad_sum_protective$sum_protective), Pheno = c(rep('Centenarian', nrow(chc_sum_protective)), rep('Offspring', nrow(off_sum_protective)), rep('Controls', nrow(ctr_sum_protective)), rep('AD', nrow(ad_sum_protective))))
    df$Pheno = factor(df$Pheno, levels = c('Centenarian', 'Offspring', 'Controls', 'AD'))
    ggplot(df, aes(x = Pheno, y = Value, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.25) + xlab('Phenotype') + ylab('Number of Protective alleles')
    df_vars = data.frame(Value = c(chc_sum_protective$n_variants_with_protective, off_sum_protective$n_variants_with_protective, ctr_sum_protective$n_variants_with_protective, ad_sum_protective$n_variants_with_protective), Pheno = c(rep('Centenarian', nrow(chc_sum_protective)), rep('Offspring', nrow(off_sum_protective)), rep('Controls', nrow(ctr_sum_protective)), rep('AD', nrow(ad_sum_protective))))
    df_vars$Pheno = factor(df_vars$Pheno, levels = c('Centenarian', 'Offspring', 'Controls', 'AD'))
    ggplot(df_vars, aes(x = Pheno, y = Value, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.25) + xlab('Phenotype') + ylab('Number of variants with at least 1 protective allele')
    df_vars_homo = data.frame(Value = c(chc_sum_protective$n_homo_protective, off_sum_protective$n_homo_protective, ctr_sum_protective$n_homo_protective, ad_sum_protective$n_homo_protective), Pheno = c(rep('Centenarian', nrow(chc_sum_protective)), rep('Offspring', nrow(off_sum_protective)), rep('Controls', nrow(ctr_sum_protective)), rep('AD', nrow(ad_sum_protective))))
    df_vars_homo$Pheno = factor(df_vars_homo$Pheno, levels = c('Centenarian', 'Offspring', 'Controls', 'AD'))
    ggplot(df_vars_homo, aes(x = Pheno, y = Value, fill = Pheno)) + geom_violin() + geom_boxplot(width=0.25) + xlab('Phenotype') + ylab('Number of variants with 2 protective allele')
    # minor allele frequency
    all_freqs
    # add protective alleles
    all_freqs$locus = paste0('chr', all_freqs$locus)
    all_freqs = merge(all_freqs, ad_snps_info, by = 'locus')
    # calculate frequency of protective alleles
    all_freqs$freq_prot_chc = ifelse(all_freqs$prot_allele == all_freqs$alt, all_freqs$freq_chc, 1-all_freqs$freq_chc)
    all_freqs$freq_prot_ctr = ifelse(all_freqs$prot_allele == all_freqs$alt, all_freqs$freq_ctr, 1-all_freqs$freq_ctr)
    all_freqs$freq_prot_ad = ifelse(all_freqs$prot_allele == all_freqs$alt, all_freqs$freq_ad, 1-all_freqs$freq_ad)
    all_freqs$freq_prot_off = ifelse(all_freqs$prot_allele == all_freqs$alt, all_freqs$freq_off, 1-all_freqs$freq_off)
    # ratio of the frequency of protective alleles
    all_freqs$ratio_freq_prot_chc_ctr = all_freqs$freq_prot_chc / all_freqs$freq_prot_ctr
    all_freqs$ratio_freq_prot_off_ctr = all_freqs$freq_prot_off / all_freqs$freq_prot_ctr
    summary(all_freqs$ratio_freq_prot_chc_ctr)
    summary(all_freqs$ratio_freq_prot_off_ctr)
    # example for apoe
    summary(chc_sum_protective[, c('chr17:44352876:C:T_C')])
    summary(ctr_sum_protective[, c('chr17:44352876:C:T_C')])
    summary(off_sum_protective[, c('chr17:44352876:C:T_C')])
    summary(ad_sum_protective[, c('chr17:44352876:C:T_C')])
    table(chc_sum_protective$'chr17:44352876:C:T_C' > 0.5)/nrow(chc_sum_protective)
    table(off_sum_protective$'chr17:44352876:C:T_C' > 0.5)/nrow(off_sum_protective)
    table(ctr_sum_protective$'chr17:44352876:C:T_C' > 0.5)/nrow(ctr_sum_protective)
    table(ad_sum_protective$'chr17:44352876:C:T_C' > 0.5)/nrow(ad_sum_protective)
    all_freqs[which(all_freqs$locus == 'chr17:44352876'),]
    table(chc_sum_protective$'chr17:44352876:C:T_C' > 1.5)/nrow(chc_sum_protective)
    table(off_sum_protective$'chr17:44352876:C:T_C' > 1.5)/nrow(off_sum_protective)
    table(ctr_sum_protective$'chr17:44352876:C:T_C' > 1.5)/nrow(ctr_sum_protective)
    table(ad_sum_protective$'chr17:44352876:C:T_C' > 1.5)/nrow(ad_sum_protective)

# PLOTS
# 1. figure 1 should be the maf, effect-size and ratio plots
    figure1_custom(singleAssoc, ratios, ad_snps, all_freqs_plot, snps_info)

# 2. figure 2 should be the density plot
    figure_2(prs_allSNPs, prs_allSNPs_noAPOE, assoc_all_combined, pheno)

# 3. figure 3 should be the simulation analysis
    figure_3(res_all_combined_nosampl, ad_snps)

# 4. figure 4 is the gene-set enrichment results -- final image is manually assembled separately
    figure_4c()

# Supplementary figures
# 1. figure s1 - raw effect-size of AD snps (AD vs. centenarians)
    figure_s1 = function(ratios, ad_snps){
    ad_cent = ratios[which(ratios$test == 'ad_vs_chc'),]
    ad_snps$locus = stringr::str_replace_all(ad_snps$locus, 'chr', '')
    ad_cent_annot = merge(ad_cent, ad_snps, by = 'locus')
    # duplicated genes
    ad_cent_annot$gene_name = ad_cent_annot$gene.x
    ad_cent_annot$gene_name[duplicated(ad_cent_annot$gene_name)] = paste0(ad_cent_annot$gene_name[duplicated(ad_cent_annot$gene_name)], ' (2)')
    # sort by ratio
    ad_cent_annot = ad_cent_annot[order(ad_cent_annot$ratio_beta),]
    ad_cent_annot$gene_name = factor(ad_cent_annot$gene_name, levels = ad_cent_annot$gene_name)
    # reshape for plot
    newdf = data.frame()
    for (i in 1:nrow(ad_cent_annot)){
        newdf = rbind(newdf, data.frame(gene = ad_cent_annot$gene_name[i], Effect = c(ad_cent_annot$beta_aligned[i], ad_cent_annot$gwas_beta[i]), Study = c('This Study', 'Bellenguez et al.')))
    }  
    # ggplot
    pdf('Downloads/figure_s1_centenarians.pdf', height = 7, width = 12)
    ggplot(data = newdf, aes(x = gene, y = Effect, fill = Study)) + geom_bar(position="dodge", stat="identity") + 
        ylim(-2.5, 2.5) + theme(legend.position = 'top', axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10)) + xlab('')
    dev.off()
    }

# 2. figure s2 - raw effect size of AD snps (AD vs. controls) + ratio
    figure_s2 = function(ratios, ad_snps){
    ad_ctr = ratios[which(ratios$test == 'ad_vs_ctr'),]
    ad_snps$locus = stringr::str_replace_all(ad_snps$locus, 'chr', '')
    ad_ctr_annot = merge(ad_ctr, ad_snps, by = 'locus')
    # duplicated genes
    ad_ctr_annot$gene_name = ad_cent_annot$gene.x
    ad_ctr_annot$gene_name[duplicated(ad_ctr_annot$gene_name)] = paste0(ad_ctr_annot$gene_name[duplicated(ad_ctr_annot$gene_name)], ' (2)')
    # sort by ratio
    ad_ctr_annot = ad_ctr_annot[order(ad_ctr_annot$ratio_beta),]
    ad_ctr_annot$gene_name = factor(ad_ctr_annot$gene_name, levels = ad_ctr_annot$gene_name)
    # reshape for plot
    newdf = data.frame()
    for (i in 1:nrow(ad_ctr_annot)){
        newdf = rbind(newdf, data.frame(gene = ad_ctr_annot$gene_name[i], Effect = c(ad_ctr_annot$beta_aligned[i], ad_ctr_annot$gwas_beta[i]), Study = c('This Study', 'Bellenguez et al.')))
    }  
    # ggplot
    fig1 = ggplot(data = newdf, aes(x = gene, y = Effect, fill = Study)) + geom_bar(position="dodge", stat="identity") + 
        ylim(-1.5, 1.5) + theme(legend.position = 'top', axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 10)) + xlab('') + ggtitle('AD cases vs. Age-matched controls')
    fig2 = ggplot(data = ad_ctr_annot, aes(x = gene_name, y = ratio_beta)) + geom_bar(stat = 'identity') + 
        xlab('') + ylab('Effect-size Ratio') + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10)) + geom_hline(yintercept = 1, linetype = 'dashed', col = 'red')
    fig3 = ggarrange(plotlist = list(fig1, fig2), nrow = 2, ncol = 1, labels = 'AUTO')
    pdf('Downloads/figure_s2_centenarians.pdf', height = 7, width = 12)
    fig3
    dev.off()
    }

# 3. figure s3 - raw effect size of AD snps (AD vs. controls) + ratio
    figure_s3 = function(ratios, ad_snps){
    ctr_chc = ratios[which(ratios$test == 'ctr_vs_chc'),]
    ad_snps$locus = stringr::str_replace_all(ad_snps$locus, 'chr', '')
    ctr_chc_annot = merge(ctr_chc, ad_snps, by = 'locus')
    # duplicated genes
    ctr_chc_annot$gene_name = ctr_chc_annot$gene.x
    ctr_chc_annot$gene_name[duplicated(ctr_chc_annot$gene_name)] = paste0(ctr_chc_annot$gene_name[duplicated(ctr_chc_annot$gene_name)], ' (2)')
    # sort by ratio
    ctr_chc_annot = ctr_chc_annot[order(ctr_chc_annot$ratio_beta),]
    ctr_chc_annot$gene_name = factor(ctr_chc_annot$gene_name, levels = ctr_chc_annot$gene_name)
    # reshape for plot
    newdf = data.frame()
    for (i in 1:nrow(ctr_chc_annot)){
        newdf = rbind(newdf, data.frame(gene = ctr_chc_annot$gene_name[i], Effect = c(ctr_chc_annot$beta_aligned[i], ctr_chc_annot$gwas_beta[i]), Study = c('This Study', 'Bellenguez et al.')))
    }  
    # ggplot
    fig1 = ggplot(data = newdf, aes(x = gene, y = Effect, fill = Study)) + geom_bar(position="dodge", stat="identity") + 
        ylim(-1.5, 1.5) + theme(legend.position = 'top', axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size = 10)) + xlab('') + ggtitle('Age-matched controls vs. Cognitively Healthy Centenarians')
    fig2 = ggplot(data = ctr_chc_annot, aes(x = gene_name, y = ratio_beta)) + geom_bar(stat = 'identity') + 
        xlab('') + ylab('Effect-size Ratio') + ylim(-6,6) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10)) + geom_hline(yintercept = 1, linetype = 'dashed', col = 'red')
    fig3 = ggarrange(plotlist = list(fig1, fig2), nrow = 2, ncol = 1, labels = 'AUTO')
    pdf('Downloads/figure_s3_centenarians.pdf', height = 7, width = 12)
    fig3
    dev.off()
    }

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

# 4. make table S5 -- same as above but comparing ad and controls
    ratios_table_s5 = ratios[which(ratios$test == 'ad_vs_ctr'),]
    ratios_table_s5$gene[which(duplicated(ratios_table_s5$gene) == TRUE)] = paste0(ratios_table_s5$gene[which(duplicated(ratios_table_s5$gene) == TRUE)], ' (2)')
    ratios_table_s5$p_fdr = p.adjust(ratios_table_s5$p, 'fdr')
    ad_snps_proxy$locus = paste(ad_snps_proxy$chrom, ad_snps_proxy$pos, sep = ":"); snps_info = ad_snps_proxy[, c('locus', 'rsid', 'maf')]
    ratios_annot_s5 = merge(ratios_table_s5, snps_info, by = 'locus')
    # reorder columns for excel
    ratios_annot_s5_final = ratios_annot_s5[, c('rsid', 'locus', 'chrom', 'pos', 'maf', 'gene', 'gwas_allele', 'gwas_beta', 'se_gwas', 'beta_aligned', 'se', 'p', 'ratio_beta', 'low_ci_ratio', 'up_ci_ratio', 'two_tail_p')]
    write.table(ratios_annot_s5_final, 'table_s5.txt', quote=F, row.names=F, sep="\t", dec=',')

# 5. make table S6 -- same as above but comparing controls and centenarians
    ratios_table_s6 = ratios[which(ratios$test == 'ctr_vs_chc'),]
    ratios_table_s6$gene[which(duplicated(ratios_table_s6$gene) == TRUE)] = paste0(ratios_table_s6$gene[which(duplicated(ratios_table_s6$gene) == TRUE)], ' (2)')
    ratios_table_s6$p_fdr = p.adjust(ratios_table_s6$p, 'fdr')
    ad_snps_proxy$locus = paste(ad_snps_proxy$chrom, ad_snps_proxy$pos, sep = ":"); snps_info = ad_snps_proxy[, c('locus', 'rsid', 'maf')]
    ratios_annot_s6 = merge(ratios_table_s6, snps_info, by = 'locus')
    # reorder columns for excel
    ratios_annot_s6_final = ratios_annot_s6[, c('rsid', 'locus', 'chrom', 'pos', 'maf', 'gene', 'gwas_allele', 'gwas_beta', 'se_gwas', 'beta_aligned', 'se', 'p', 'ratio_beta', 'low_ci_ratio', 'up_ci_ratio', 'two_tail_p')]
    write.table(ratios_annot_s6_final, 'table_s.txt', quote=F, row.names=F, sep="\t", dec=',')

# 6. make table s7 -- prs associations
    write.table(assoc_all_combined, 'table_s7.txt', quote=F, row.names=F, sep='\t', dec=',')

# WORKSPACE
    save.image('path/to/workspace.RData')

##################################################################
# REVIEWER'S COMMENTS
    # Synergistic effect of SNPs
        # load environment
            load('path/to/workspace.RData')
        # combination of 86 snps in groups of 2
            snplist_all = colnames(dosages_pheno)[grep(':', colnames(dosages_pheno))]
            combinations_snps <- combn(snplist_all, 2)
            list_of_vectors <- lapply(1:ncol(combinations_snps), function(i) combinations_snps[, i])
            # iterate through snps
            interaction_models = data.frame()
            for (pair in list_of_vectors){
                # make subset of the dosage file
                sb = dosages_pheno[, pair]; sb$PHENO = dosages_pheno$diagnosis; sb$ID = dosages_pheno$IID;
                colnames(sb) = c('snp1', 'snp2', 'pheno', 'id')
                sb = merge(sb, pheno, by.x = 'id', by.y = 'ID_GWAS')
                # add case-control labels
                sb_ad = sb[which(sb$pheno %in% c("Centenarian", "Probable_AD", "Possible_AD", "AD_path")),]
                sb_ctr = sb[which(sb$pheno %in% c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
                sb_ad$case_control_status = ifelse(sb_ad$pheno == 'Centenarian', 0, 1)
                sb_ctr$case_control_status = ifelse(sb_ctr$pheno == 'Centenarian', 0, 1)
                # make the model
                model_ad = glm(case_control_status ~ snp1 * snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ad, family = 'binomial')
                model_ctr = glm(case_control_status ~ snp1 * snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ctr, family = 'binomial')
                # save info
                if ('snp1:snp2' %in% rownames(summary(model_ad)$coefficients)){
                    interaction_models = rbind(interaction_models, data.frame(snp_pair = paste(pair, collapse="-"), beta_interaction = summary(model_ad)$coefficients['snp1:snp2', c("Estimate")], se_interaction = summary(model_ad)$coefficients['snp1:snp2', c('Std. Error')], p_interaction = summary(model_ad)$coefficients['snp1:snp2', c('Pr(>|z|)')], pheno = 'AD'))
                    interaction_models = rbind(interaction_models, data.frame(snp_pair = paste(pair, collapse="-"), beta_interaction = summary(model_ctr)$coefficients['snp1:snp2', c("Estimate")], se_interaction = summary(model_ctr)$coefficients['snp1:snp2', c('Std. Error')], p_interaction = summary(model_ctr)$coefficients['snp1:snp2', c('Pr(>|z|)')], pheno = 'CTR'))
                }
            }
            # separate ad from ctr associations
            interaction_models_ad = interaction_models[which(interaction_models$pheno == 'AD'),]
            interaction_models_ctr = interaction_models[which(interaction_models$pheno == 'CTR'),]
            # correct with fdr
            interaction_models_ad$p_interaction_fdr = p.adjust(interaction_models_ad$p_interaction, 'fdr')
            interaction_models_ctr$p_interaction_fdr = p.adjust(interaction_models_ctr$p_interaction, 'fdr')
            # sort by p
            interaction_models_ad = interaction_models_ad[order(interaction_models_ad$p_interaction_fdr),]
            interaction_models_ctr = interaction_models_ctr[order(interaction_models_ctr$p_interaction_fdr),]
            head(interaction_models_ad)
            # only 1 significant after correction in AD -- TREML2 (chr6:41181270:A:G_A) and IGH-2 (chr14:106665591:G:A_G)
            # 0 carriers in centenarians, 8 in the AD cases, 7 of which are homozygous carriers of the other variant
            head(interaction_models_ctr)
            # 2 significant after correction in CTR -- TREML2 (chr6:41181270:A:G_A) and COX7C (chr5:86927378:T:C_T) // TREML2 (chr6:41181270:A:G_A) and KLF16 (chr19:1854254:G:GC_G )
            # they are all common expect TREML2 which is rare (0.004) - maybe that explains?
        # maybe interesting to look at associations with the components -- maybe not
            association_components = data.frame()
            for (snp in snplist_all){
                sb = data.frame(snp1 = dosages_pheno[, snp], pheno = dosages_pheno$diagnosis, id = dosages_pheno$IID)
                sb = merge(sb, pheno, by.x = 'id', by.y = 'ID_GWAS')
                # add case-control labels
                sb_ad = sb[which(sb$pheno %in% c("Centenarian", "Probable_AD", "Possible_AD", "AD_path")),]
                sb_ctr = sb[which(sb$pheno %in% c("Centenarian", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
                sb_mat = sb[which(sb$pheno %in% c("Probable_AD", "Possible_AD", "AD_path", "Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
                sb_ad$case_control_status = ifelse(sb_ad$pheno == 'Centenarian', 0, 1)
                sb_ctr$case_control_status = ifelse(sb_ctr$pheno == 'Centenarian', 0, 1)
                sb_mat$case_control_status = ifelse(sb_mat$pheno %in% c("Probable_AD", "Possible_AD", "AD_path"), 1, 0)
                # make the model
                model_ad = glm(case_control_status ~ snp1 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ad, family = 'binomial')
                model_ctr = glm(case_control_status ~ snp1 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ctr, family = 'binomial')
                model_mat = glm(case_control_status ~ snp1 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_mat, family = 'binomial')
                # of the PCs, take the most significant
                tmp_ad = data.frame(summary(model_ad)$coefficients); tmp_ad = tmp_ad[grep('PC', rownames(tmp_ad)),]; tmp_ad = tmp_ad[order(tmp_ad$'Pr...z..'),]
                tmp_ctr = data.frame(summary(model_ctr)$coefficients); tmp_ctr = tmp_ctr[grep('PC', rownames(tmp_ctr)),]; tmp_ctr = tmp_ctr[order(tmp_ctr$'Pr...z..'),]
                tmp_mat = data.frame(summary(model_mat)$coefficients); tmp_mat = tmp_mat[grep('PC', rownames(tmp_mat)),]; tmp_mat = tmp_mat[order(tmp_mat$'Pr...z..'),]
                # save the info
                association_components = rbind(association_components, data.frame(snp = rep(snp, 3), type = c('AD_CHC', 'CTR_CHC', 'AD_CTR'), variable = c(rownames(tmp_ad)[1], rownames(tmp_ctr)[1], rownames(tmp_mat)[1]), p = c(tmp_ad$'Pr...z..'[1], tmp_ctr$'Pr...z..'[1], tmp_mat$'Pr...z..'[1])))
            }
            # separate associations
            association_components_ad = association_components[which(association_components$type == 'AD_CHC'),]
            association_components_ctr = association_components[which(association_components$type == 'CTR_CHC'),]
            association_components_mat = association_components[which(association_components$type == 'AD_CTR'),]
            # correct with fdr
            association_components_ad$p_fdr = p.adjust(association_components_ad$p, 'fdr', n=nrow(association_components_ad))
            association_components_ctr$p_fdr = p.adjust(association_components_ctr$p, 'fdr', n=nrow(association_components))
            association_components_mat$p_fdr = p.adjust(association_components_mat$p, 'fdr', n=nrow(association_components))
            # sort by p
            association_components_ad = association_components_ad[order(association_components_ad$p_fdr),]
            association_components_ctr = association_components_ctr[order(association_components_ctr$p_fdr),]
            association_components_mat = association_components_mat[order(association_components_mat$p_fdr),]
            head(association_components_ad)
            head(association_components_ctr)
            head(association_components_mat)
        # maybe we can look at survival in the centenarians
            # read survival data
            library(survival)
            surv_data = read.table('age_death_centenarians.txt', h=T, stringsAsFactors=F, sep="\t", dec=",")
            # take centenarians and merge with this data
            chc = pheno[which(pheno$diagnosis == 'Centenarian'),]
            chc = merge(chc, pheno_final_raw, by = 'ID_GWAS')
            chc_info = chc[, c('ID_GWAS', 'ID_100plus', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'sex.x')]
            colnames(chc_info) = c('ID_GWAS', 'ID_100plus', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'sex')
            prs_chc_info = merge(prs_allSNPs, chc_info, by.x = 'IID', by.y = 'ID_GWAS')
            prs_chc_info_surv = merge(prs_chc_info, surv_data, by.x = 'ID_100plus', by.y = 'studysubjectid')
            # we lose 1 sample
            # 27 are still alive --> calculate difference in time between age at entrance and now
            tmp = prs_chc_info_surv
            xx = data.frame(stringr::str_split_fixed(tmp$date_entrance_interview, '-', 3)); xx$X3 = paste0('20', xx$X3); xx$X5 = paste0(xx$X1, '-', xx$X2, '-', xx$X3)
            tmp$date_entrance_interview = xx$X5
            tmp$date_entrance_interview = as.Date(tmp$date_entrance_interview, format = '%d-%m-%Y')
            tmp$death = ifelse(is.na(tmp$age_at_death_calc), 0, 1)
            tmp$age_at_death_calc[is.na(tmp$age_at_death_calc)] = tmp$age_entrance_interview[is.na(tmp$age_at_death_calc)] + difftime(as.Date('24-10-2023', format = '%d-%m-%Y'), tmp$date_entrance_interview[is.na(tmp$age_at_death_calc)], unit="weeks")/52.25
            info_surv = tmp[, c('IID', 'ID_100plus', 'death', 'age_at_death_calc', 'sex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
            prs_noAPOE_surv = merge(prs_allSNPs_noAPOE, info_surv, by = 'IID')
            # exclude 2 NAs
            tmp_clean = tmp[!is.na(tmp$age_at_death_calc),]
            cox_model = coxph(Surv(age_at_death_calc, death) ~ PRS + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = tmp)
            cox_model_noAPOE = coxph(Surv(age_at_death_calc, death) ~ PRS + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_noAPOE_surv)
            summary(cox_model)
            summary(cox_model_noAPOE)
            # try survival in lasa as well
            surv_data = read.table('lasa_phenotypes_QCpassed.txt', h=T, stringsAsFactors=F, sep="\t")
            surv_data$surv_time <- surv_data$age_GWA - surv_data$age_baseline
            surv_data <- surv_data[, c("respnr", "death", "surv_time", "age_baseline", "sex")]
            surv_data <- surv_data[which(surv_data$death != "-1"),]
            surv_data$respnr <- as.character(surv_data$respnr)
            # add prs
            surv_data_prs = merge(prs_allSNPs, surv_data, by.y = 'respnr', by.x = 'IID')
            surv_data_prs_noAPOE = merge(prs_allSNPs_noAPOE, surv_data, by.y = 'respnr', by.x = 'IID')
            # add covariates
            covar = pheno[, c('ID_GWAS', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5')]
            surv_data_prs = merge(surv_data_prs, covar, by.x = 'IID', by.y = 'ID_GWAS')
            surv_data_prs_noAPOE = merge(surv_data_prs_noAPOE, covar, by.x = 'IID', by.y = 'ID_GWAS')
            head(surv_data_prs)
            # cox model
            cox_lasa = coxph(Surv(time = surv_data_prs$age_baseline-min(surv_data_prs$age_baseline), time2 = surv_data_prs$age_baseline+surv_data_prs$surv_time-min(surv_data_prs$age_baseline), event = surv_data_prs$death, type="counting") ~ PRS + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = surv_data_prs)
            cox_lasa_noAPOE = coxph(Surv(time = surv_data_prs_noAPOE$age_baseline-min(surv_data_prs_noAPOE$age_baseline), time2 = surv_data_prs_noAPOE$age_baseline+surv_data_prs_noAPOE$surv_time-min(surv_data_prs_noAPOE$age_baseline), event = surv_data_prs_noAPOE$death, type="counting")  ~ PRS + sex + PC1 + PC2 + PC3 + PC4 + PC5, data = surv_data_prs_noAPOE)
            summary(cox_lasa)
            summary(cox_lasa_noAPOE)
    # Compare offspring vs. controls -- single-variant
        single_variant_offspring = data.frame()
        for (snp in snplist_all){
            sb = data.frame(snp1 = dosages_pheno[, snp], pheno = dosages_pheno$diagnosis, id = dosages_pheno$IID)
            colnames(sb) = c('snp1', 'pheno', 'id')
            # add case-control labels
            sb_off = sb[which(sb$id %in% children_clean$IID),]; sb_off$case_control_status = 0
            sb_off = sb_off[which(sb_off$pheno == 'family_100plus'),]
            sb_ctr = sb[which(sb$pheno %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]; sb_ctr$case_control_status = 1
            sb = rbind(sb_off, sb_ctr)
            sb_ph = merge(sb, children_clean, by.x = 'id', by.y = 'IID')
            # make the model
            model = glm(case_control_status ~ snp1 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ph, family = 'binomial')
            # of the PCs, take the most significant
            tmp = data.frame(summary(model)$coefficients)
            single_variant_offspring = rbind(single_variant_offspring, data.frame(snp = snp, beta = tmp$Estimate[2], se = tmp$"Std..Error"[2], p = tmp$"Pr...z.."[2]))        
        }
        # correct
        single_variant_offspring$p_fdr = p.adjust(single_variant_offspring$p, 'fdr')
        # sort
        single_variant_offspring = single_variant_offspring[order(single_variant_offspring$p_fdr),]
        # add snps info
        single_variant_offspring$snpid = NA
        for (i in 1:nrow(single_variant_offspring)){ single_variant_offspring$snpid[i] = paste0(stringr::str_split_fixed(single_variant_offspring$snp[i], ':', 4)[, 1:2], collapse=":") }
        head(single_variant_offspring)
        single_variant_offspring = merge(single_variant_offspring, ad_snps, by.x = 'snpid', by.y = 'locus')
        single_variant_offspring$allele = stringr::str_split_fixed(stringr::str_split_fixed(single_variant_offspring$snp, ':', 4)[, 4], '_', 2)[, 2]
        write.table(single_variant_offspring, '20231023_single_variant_offspring.txt', quote=F, row.names=F, sep="\t", dec=',')
    # Sensitivity analysis for EOAD (971, age <=65) vs. LOAD (1309, age>65)
        single_variant_sensitivity = data.frame()
        for (snp in snplist_all){
            sb = data.frame(snp1 = dosages_pheno[, snp], pheno = dosages_pheno$diagnosis, id = dosages_pheno$IID)
            colnames(sb) = c('snp1', 'pheno', 'id')
            # only ad
            sb_ad = sb[which(sb$pheno %in% c("Probable_AD", "Possible_AD", "AD_path")),]
            # add case-control labels
            sb_ad_ph = merge(sb_ad, pheno, by.x = 'id', by.y = 'ID_GWAS')
            # add case-control status based on the age
            sb_ad_ph$case_control_status = ifelse(sb_ad_ph$age <=65, 1, 0)
            # make the model
            model = glm(case_control_status ~ snp1 + PC1 + PC2 + PC3 + PC4 + PC5, data = sb_ad_ph, family = 'binomial')
            # of the PCs, take the most significant
            tmp = data.frame(summary(model)$coefficients)
            single_variant_sensitivity = rbind(single_variant_sensitivity, data.frame(snp = snp, beta = tmp$Estimate[2], se = tmp$"Std..Error"[2], p = tmp$"Pr...z.."[2]))        
        }
        # correct
        single_variant_sensitivity$p_fdr = p.adjust(single_variant_sensitivity$p, 'fdr')
        # sort
        single_variant_sensitivity = single_variant_sensitivity[order(single_variant_sensitivity$p_fdr),]
        head(single_variant_sensitivity)
        # add snps info
        single_variant_sensitivity$snpid = NA
        for (i in 1:nrow(single_variant_sensitivity)){ single_variant_sensitivity$snpid[i] = paste0(stringr::str_split_fixed(single_variant_sensitivity$snp[i], ':', 4)[, 1:2], collapse=":") }
        head(single_variant_sensitivity)
        single_variant_sensitivity = merge(single_variant_sensitivity, ad_snps, by.x = 'snpid', by.y = 'locus')
        single_variant_sensitivity$allele = stringr::str_split_fixed(stringr::str_split_fixed(single_variant_sensitivity$snp, ':', 4)[, 4], '_', 2)[, 2]
        write.table(single_variant_sensitivity, '20231023_single_variant_sensitivity_analysis.txt', quote=F, row.names=F, sep="\t", dec=',')
        # and also at the PRS level
        prs_ad = prs_allSNPs[which(prs_allSNPs$PHENO %in% c("Probable_AD", "Possible_AD", "AD_path")),]
        prs_ad_noAPOE = prs_allSNPs_noAPOE[which(prs_allSNPs_noAPOE$PHENO %in% c("Probable_AD", "Possible_AD", "AD_path")),]
        prs_ad = merge(prs_ad, pheno, by.x = "IID", by.y = 'ID_GWAS')
        prs_ad_noAPOE = merge(prs_ad_noAPOE, pheno, by.x = "IID", by.y = 'ID_GWAS')
        # separate eoad and load
        prs_ad$case_control_status = ifelse(prs_ad$age <= 65, 1, 0)
        prs_ad_noAPOE$case_control_status = ifelse(prs_ad_noAPOE$age <= 65, 1, 0)
        # models
        sensitivity_model = glm(case_control_status ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_ad)
        sensitivity_model_noAPOE = glm(case_control_status ~ PRS + PC1 + PC2 + PC3 + PC4 + PC5, data = prs_ad_noAPOE)
        summary(sensitivity_model)
        summary(sensitivity_model_noAPOE)
    # Additional cohort info and plot
        ad = pheno[which(pheno$diagnosis %in% c('Probable_AD', 'Possible_AD', 'AD_path')),]
        ctr = pheno[which(pheno$diagnosis %in% c("Control_100plus", "Control_LASA", "Control_other_twin", "Control_path", "SCD")),]
        chc = pheno[which(pheno$diagnosis %in% c('Centenarian')),]
        off = children_clean[which(children_clean$PHENO %in% c('family_100plus')),]
        summary(ad$age)
        summary(ctr$age)
        summary(chc$age)
        summary(pheno_final_raw$age[which(pheno_final_raw$ID_GWAS %in% off$IID)])
        off2 = pheno_final_raw[which(pheno_final_raw$ID_GWAS %in% off$IID),]
        library(ggplot2)
        data_ages = data.frame(Type = c(rep('AD', nrow(ad)), rep('Controls', nrow(ctr)), rep('Centenarians', nrow(chc))), age = c(ad$age, ctr$age, chc$age))
        # rename classes
        data_ages$Type[which(data_ages$Type == 'AD')] = 'AD cases'
        data_ages$Type[which(data_ages$Type == 'Controls')] = 'Age-matched controls'
        data_ages$Type[which(data_ages$Type == 'Centenarians')] = 'Cognitively Healthy Centenarians'
        data_ages$Type = factor(data_ages$Type, levels = c('AD cases', 'Age-matched controls', 'Cognitively Healthy Centenarians'))
        plt = ggplot(data_ages, aes(x = Type, y = age, fill = Type)) + geom_violin(alpha = 0.5, width = 1.4) + geom_boxplot(width = 0.07) + theme_bw() + theme(legend.position = 'top', axis.text = element_text(size = 16), axis.title = element_text(size = 18), legend.text = element_text(size = 14)) + xlab('Sample type') + ylab('Age at onset/inclusion')
        pdf('Downloads/figure_s1_ages.pdf', height=8, width = 14)
        print(plt)
        dev.off()
    # Save workspace updated
        save.image('path/to/updated_workspace.RData')
        
        
        