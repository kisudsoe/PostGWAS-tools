#!/usr/bin/env Rscript

## Help Messages ##
help_message = '
Version: 2020-08-14

Usgae:
    Rscript ./src/venn_analysis.r --table --base <Comparing BED file paths> --out <out folder path> [option]
    Rscript ./src/venn_analysis.r --fig --base <Comparing BED file paths> --out <out folder path> [option]

Functions:
    --table  Generating venn table
    --fig    Generating venn diagram plot

Opions:
    --title  Venn plot title
'


## Command Input ##
commandLine = commandArgs(trailingOnly=T)
comm = paste(unlist(commandLine),collapse=' ')
listoptions = unlist(strsplit(comm,'--'))[-1]
args = lapply(listoptions,function(x) {
    arg = unlist(strsplit(x,' '))[-1]
    if(length(arg)==0) { return(TRUE)
    } else return(arg)
})
args_names = sapply(listoptions,function(x) {
    option = unlist(strsplit(x,' '))[1]
})
args_names = unlist(args_names)
names(args) = args_names


## Function Start ##
venn_tb = function(
    grouplist = NULL,
    out       = NULL,
    title     = NULL
) {
    wfig  = FALSE
    wfile = TRUE
    venn_analysis(
        grouplist = grouplist,
        out       = out,
        title     = title,
        wfig      = wfig,
        wfile     = wfile,
        verbose   = FALSE
    )
}


venn_fig = function(
    grouplist = NULL,
    out       = NULL,
    title     = NULL
) {
    wfig  = TRUE
    wfile = FALSE
    venn_analysis(
        grouplist = grouplist,
        out       = out,
        title     = title,
        wfig      = wfig,
        wfile     = wfile,
        verbose   = FALSE
    )
}


venn_analysis = function(
    grouplist = NULL,
    out       = NULL,
    title     = NULL,
    wfig      = TRUE,
    wfile     = TRUE,
    verbose   = FALSE
) {
    g_names = names(grouplist)
    unionlist = Reduce(union,grouplist)
    
    # Collapse the list to one dataFrame list
    unionPr=NULL; g_titles=NULL
    for(i in 1:length(grouplist)) {
        unionPr = cbind(unionPr,grouplist[[i]][match(unionlist,grouplist[[i]])])
        g_titles = c(g_titles,paste0(g_names[i],'\n',length(grouplist[[i]])))
    }
    rownames(unionPr) = unionlist
    
    # Generate binary table to match with ID list
    union = (unionPr != '')			# Transcform values to TRUE, if ID exists.
    union[is.na(union)] = FALSE		# Transform NA to FALSE value
    union = as.data.frame(union)	# Make 'union' to data.frame from
    colnames(union) = g_titles		# Names attach to venn diagram
    
    # Draw Venn Diagram
    if(is.null(title)) title = paste0('Venn analysis of ',nrow(union),' genes')
    if(wfig & length(grouplist)==3) {
        suppressMessages(library(eulerr))
        suppressMessages(library(limma))
        f.name1 = paste0(out,'/venn_',g_names[1],'_',g_names[2],'.png')
        venn_counts = vennCounts(union)
        venn_c = unlist(venn_counts[,4])
        png(f.name1,width=10,height=10,units='in',res=100)
        # ref: https://creately.com/blog/diagrams/venn-diagrams-vs-euler-diagrams/
        venn_fit = euler(c(
            'A'    =as.numeric(venn_c[5]), # How to automate this?!
            'B'    =as.numeric(venn_c[3]),
            'C'    =as.numeric(venn_c[2]),
            'A&B'  =as.numeric(venn_c[7]),
            'A&C'  =as.numeric(venn_c[6]),
            'B&C'  =as.numeric(venn_c[4]),
            'A&B&C'=as.numeric(venn_c[8])
                    )) # eulerr
        print(venn_fit)
        cat(paste0('\nEuler fit is done.\n'))
        p=plot(
            venn_fit, quantities=T,
            labels=colnames(venn_counts)[1:3],
            edges=list(col=c('red','green','blue'),lwd=3),
            fills=list(fill=c('red','green','blue'),alpha=0.2),
            main=''
            ) # eulerr
        print(p)
        graphics.off() # killing all devices
        cat(paste0('\nFigure draw: ',f.name1,'\n'))
    } else if(wfig & length(grouplist)%in%c(2,4:6)) {
        suppressMessages(library(limma))
        f.name1 = paste0(out,'/venn_',g_names[1],'_',g_names[2],'.png')
        png(f.name1,width=10,height=10,units='in',res=100)
        vennDiagram(union, main=title, circle.col=rainbow(length(g_titles))) # limma
        graphics.off() # killing all devices
        cat(paste0('\nFigure draw: ',f.name1,'\n'))
    } else {
        if(verbose) paste0('[NOTICE] Too many input groups..\n') %>% cat
    }

    # Write files
    colnames(union) = g_names
    union_df = cbind(union,ids=rownames(union))
    if(wfile) {
        f.name2  = paste0(out,'/venn_',g_names[1],'_',g_names[2],'.tsv')
        write.table(union_df,f.name2,row.names=F,col.names=T,quote=F,sep='\t')
        cat(paste0('File write: ',f.name2,'\n'))
    }
    return()
}


base2list = function(
    b_paths = NULL
) {
    suppressMessages(library(dplyr))
    suppressMessages(library(data.table))
    suppressMessages(library(tools))

    out_li = lapply(b_paths,function(x) {
        x_ext = file_ext(x)
        if(x_ext=='bed') {
            df = fread(x,header=F)[,4] %>% unlist
        } else {
            paste0('[Error] ',x,' is not BED file.\n') %>% cat
            return(NULL)
        }
    })
    base_nms = tools::file_path_sans_ext(b_paths %>% basename)
    names(out_li) = base_nms
    return(out_li)
}


## Parsing Arguments ##
if(length(args$help)>0) {  help    = args$help
} else                     help    = FALSE
if(help) { cat(help_message); quite() }

# Functions
if(length(args$table)>0) { opt_tb  = args$table
} else                     opt_tb  = FALSE
if(length(args$fig)>0) {   opt_fig = args$fig
} else                     opt_fig = FALSE

# Arguments
if(length(args$base)>0)    b_path  = args$base
if(length(args$out)>0)     out     = args$out
if(length(args$title)>0) { title   = args$title
} else                     title   = NULL

# Get list
b_list = base2list(b_path)
if(is.null(b_list)) {
    paste0('\n[Error] NULL list returned. Nothing read.\n') %>% cat
    return()
}

## Run function
if(opt_tb) {
    venn_tb(b_list,out,title)
} else if(opt_fig) {
    venn_fig(b_list,out,title)
}