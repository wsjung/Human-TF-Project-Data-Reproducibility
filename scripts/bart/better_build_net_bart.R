generate_bart_net = function(p_in_expr_target
                             , p_in_expr_reg
                             , p_out_dir
                             , ntree
                             , seed
                             , p_src_code
			     , batch_file
                            ){
    
    df_expr_target = read.csv(p_in_expr_target, header=TRUE, row.names=1, sep="\t", numerals = "no.loss")
    print(tail(row.names(df_expr_target)))
    batch_genes = scan(batch_file, what="", sep="\n")
    print(tail(batch_genes))
    df_expr_target = df_expr_target[batch_genes,] # SUBSET DF_EXPR_TARGET TO ONLY INCLUDE BATCH GENES
    gc() # RUN GARBAGE COLLECTOR TO SAVE MEMORY
    print(tail(row.names(df_expr_target)))

    row_names <- row.names(df_expr_target)
    df_expr_reg = read.csv(p_in_expr_reg, header=TRUE, row.names=1, sep="\t", numerals = "no.loss")
    l_in_target = as.factor(rownames(df_expr_target))
    l_in_reg = as.factor(rownames(df_expr_reg))
    l_in_sample = as.factor(colnames(df_expr_target))

    options(digits=22)
    df_expr_reg = data.frame(lapply(df_expr_reg, function(x) as.double(as.character(x))))
    df_expr_target = data.frame(lapply(df_expr_target, function(x) as.double(as.character(x))))
    
    # generate intermediate files (allowed and perturbed)
    df_allowed_perturbed = generate_allowed_perturbed_matrices(l_in_target, l_in_reg, l_in_sample, NULL, p_src_code)
    df_allowed = sapply(as.data.frame(as.matrix(df_allowed_perturbed[[1]])), as.logical)
    df_perturbed = sapply(as.data.frame(as.matrix(df_allowed_perturbed[[2]])), as.logical)
    
    # masking pertubed entries in response from training
    df_expr_target[df_perturbed] = NA
    row.names(df_expr_target) <- row_names
    print(tail(row.names(df_expr_target)))
    df_bart_net = getBartNetwork(tgtLevel=t(as.matrix(df_expr_target))
				 , tfLevel=t(as.matrix(df_expr_reg))
                                 , regMat=df_allowed
                                 , ntree=ntree)
    					
    write.table(df_bart_net$regScore
                , file.path(p_out_dir, basename(batch_file))
                , row.names=l_in_reg
                , col.names=l_in_target
                , quote=FALSE
                , sep="\t")
}

if (sys.nframe() == 0){
    # =========================================== #
    # |       *** Install packages ***          | #
    # =========================================== #
    library("optparse")
    
    # =========================================== #
    # |         **** Parse Arguments ****       | #
    # =========================================== #
    p_in_expr_target = make_option(c("--p_in_expr_target"), type="character", help='input - path of expression of target genes', default=NULL)
    p_in_expr_reg = make_option(c("--p_in_expr_reg"), type="character", help="input - path of expression of regulators", default=NULL)
    p_out_dir = make_option(c("--p_out_dir"), type="character", default=NULL, help="output - path of output directory for results")
    ntree = make_option(c("--ntree"), type="integer", help="number of trees for BART algo")
    seed = make_option(c("--seed"), type="integer", default=747, help="seed for reproducibility")
    p_src_code = make_option(c("--p_src_code"), type="character", default=NULL, help="path of the source code")
    p_scripts = make_option(c("--p_scripts"), type="character", default=NULL, help="path of Daniel's updated scripts")
    batch_file = make_option(c("--batch_file"), type="character", help="batch of target genes to analyze")
    opt_parser = OptionParser(option_list=list(p_in_expr_target, p_in_expr_reg, p_out_dir, ntree, seed, p_src_code, p_scripts, batch_file))
    opt = parse_args(opt_parser)
    
    if (is.null(opt$p_in_expr_target) || is.null(opt$p_in_expr_reg) || is.null(opt$p_out_dir) || is.null(opt$ntree) || is.null(opt$p_src_code) || is.null(opt$p_scripts) || is.null(opt$batch_file)
       )
    {
        print_help(opt_parser)
        stop("Arguments p_in_expr_target, p_in_expr_reg, p_out_dir, ntree, p_src_code, p_scripts, batch_file are mandatory")
    }
    
    # load local R
    source(paste(opt$p_scripts, "better_build_bart_network.r", sep=""))
    source(paste(opt$p_src_code, "src/build_lasso/code/prepare_data_generate_allowed_perturbed_and_scale_normalize.R", sep=""))
    
    quit(status=generate_bart_net(p_in_expr_target=opt$p_in_expr_target
                                  , p_in_expr_reg=opt$p_in_expr_reg
                                  , p_out_dir=opt$p_out_dir
                                  , ntree=opt$ntree
                                  , seed=opt$seed
                                  , p_src_code=opt$p_src_code
				  , batch_file=opt$batch_file
                                 ))
}
