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
        cat(paste0('> Euler fit is done.'))
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
    return(union_df)
}