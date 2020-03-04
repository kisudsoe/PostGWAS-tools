help_message = '
db_venn, v2020-03-01
This is a function call for venn analysis of filtered DB data.

Usage: Rscript postgwas-exe.r --dbvenn <function> --base <base files> --out <out folder> --fig <figure out folder>

Functions:
	venn    Venn analysis of rsids.
	summ    Generating summary table with annotations.

Global arguments:
	--base  <base files>
			At least 2 Base BED file paths are mendatory.
			If the base file number is over 4, then venn diagram is not generated.
	--out   <out folder>
			Out folder path is mendatory. Default is "data" folder.

Required arguments:
	--fig   <figure out folder>
			An optional argument for the "venn" function to save figure file.
			If no figure out path is designated, no venn figure will generated.
	--uni   <defulat: FALSE>
			An optional argument for the "venn" function to save union BED file.
	--ann   <annotation files>
			An optional argument for the "summ" function. Add annotations to the summary table.
'

## Load global libraries ##
suppressMessages(library(dplyr))

## Functions Start ##
summ_ann = function(
	f_paths   = NULL,  # Input BED file paths
	ann_paths = NULL,  # Add annotation file paths
	out       = 'data' # Out folder path
) {
	# Prepare...
	paste0('\n** Run function: db_venn.r/venn...\n') %>% cat
	quit()
}

venn_bed = function(
	f_paths  = NULL,   # Input BED file paths
	out      = 'data', # Out folder path
	fig      = NULL,   # Figure out folder path
	uni_list = FALSE   # 
) {
	# Function specific library
	suppressMessages(library(eulerr))
	suppressMessages(library(tools))
	
	# Prepare...
	paste0('\n** Run function: db_venn.r/venn_bed...\n') %>% cat
	n         = length(f_paths)
	snp_li    = list()
	snpids_li = list()
	names     = NULL

	# Read BED files
	for(i in 1:n) {
		tb = read.delim(f_paths[i],header=F)
		colnames(tb)   = c('chr','start','end','ann')
		snp_li[[i]]    = tb
		snpids_li[[i]] = tb$ann
		f_base = basename(f_paths[i]) %>% file_path_sans_ext
		names = c(names,f_base)
		paste0('  Read: ',f_base,'\n') %>% cat
	}
	snp_df    = data.table::rbindlist(snp_li) %>% unique
	unionlist = Reduce(union,snpids_li)

	# Prepare for collapsing the list to one dataFrame list
	unionPr  = NULL
	subtitle = NULL
	for(i in 1:length(snpids_li)) {
		unionPr  = cbind(unionPr,snpids_li[[i]][match(unionlist,snpids_li[[i]])])
		subtitle = c(subtitle,paste0(names[i],'\n',length(snpids_li[[i]])))
	}
	rownames(unionPr) = unionlist

	# Generate binary table to match with ID list
	union = (unionPr != '')			# Transcform values to TRUE, if ID exists.
	union[is.na(union)] = FALSE		# Transform NA to FALSE value
	union = as.data.frame(union)	# Make 'union' to data.frame from
	colnames(union) = subtitle		# Names attach to venn diagram

	# Draw Venn plot
	title  = paste0(' Venn analysis of ',nrow(union),' SNPs')
	if(n %in% c(2,4) & !is.null(fig)) {
		fig_name1 = paste0(fig,'/venn_',n,'_',paste0(names,collapse=', '),'.png')
		png(fig_name1,width=10,height=10,units='in',res=100)
		limma::vennDiagram(
			union,
			main       = title,
			circle.col = rainbow(subtitle%>%length)
		)
		dev.off()
		paste0('\nFigure draw:\t\t',fig_name1,'\n') %>% cat
	} else if(n>4) message("\n[Message] Can't plot Venn diagram for more than 5 sets.")

	# Draw Euler plot
	if(n == 3 & !is.null(fig)) {
		paste0('\n** Euler fitting... ') %>% cat
		vennCount = limma::vennCounts(union)
		fig_name2 = paste0(fig,'/euler_',n,'_',paste0(names,collapse=', '),'.png')
		png(fig_name2,width=10,height=10,units='in',res=100)
		venn_c   = unlist(vennCount[,4])
		venn_fit = euler(c(
			'A'     = venn_c[5] %>% as.numeric,
			'B'     = venn_c[3] %>% as.numeric,
			'C'     = venn_c[2] %>% as.numeric,
			'A&B'   = venn_c[7] %>% as.numeric,
			'A&C'   = venn_c[6] %>% as.numeric,
			'B&C'   = venn_c[4] %>% as.numeric,
			'A&B&C' = venn_c[8] %>% as.numeric
		))
		paste0('done.\n') %>% cat
		p = plot(
			venn_fit,
			quantities = T,
			labels     = colnames(vennCount)[1:3],
			edges      = list(col  = c('red','green','blue'),lwd=3),
			fills      = list(fill = c('red','green','blue'),alpha=.2),
			main       = title
		)
		print(p)
		dev.off()
		paste0('\nFigure draw:\t\t',fig_name2,'\n') %>% cat
	} else message("\n[Message] Can't plot Euler plot.")

	# Save the venn result as a TSV file
	colnames(union) = names
	union_df  = cbind(union,ann=rownames(union))
	union_out = merge(snp_df,union_df,by='ann',all=T)
	if(n>3) {
		con_name = paste0(names[1],'-',names[n])
	} else con_name = paste0(names,collapse=', ')
	f_name1   = paste0(out,'/venn_',n,'_',con_name,'.tsv')
	write.table(union_out,f_name1,row.name=F,quote=F,sep='\t')
	paste0('Write TSV file:\t\t',f_name1,'\n') %>% cat

	# Write core snp list as a BED file
	if(n == 2 & !uni_list) {
		which_row = which(
			union_out[,5] == TRUE &
			union_out[,6] == TRUE
		)
	} else if(n == 3 & !uni_list) {
		which_row = which(
			union_out[,5] == TRUE &
			union_out[,6] == TRUE &
			union_out[,7] == TRUE
		)
	} else if(uni_list) {
		which_row = c(1:nrow(union_out)) # all union list
	} else {
		paste0('\n[Message] If you need a core rsid BED file,
	please input two or three files.') %>% message
		quit()
	}
	core_df = union_out[which_row, c(2:4,1)]
	f_name2 = paste0(out,'/snp_core_',unique(core_df)%>%nrow,'.bed')
	write.table(core_df,f_name2,row.names=F,col.names=F,quote=F,sep='\t')
	paste0('Write core BED file:\t',f_name2,'\n\n') %>% cat
}

db_venn = function(
	args = NULL
) {
	# Get help
	if(length(args$help)>0) {   help     = args$help
    } else                      help     = FALSE
    if(help) {                  cat(help_message); quit() }

	# Global arguments
	if(length(args$base)>0)     b_path   = args$base
    if(length(args$out)>0)      out      = args$out
    if(length(args$debug)>0) {  debug    = args$debug
    } else                      debug    = FALSE

	# Reguired arguments
	if(length(args$fig)>0)      fig_path = args$fig
	if(length(args$uni)>0) {    uni_list = args$uni
	}                           uni_list = FALSE
	if(length(args$ann)>0)      anns     = args$ann

	# Run function
	source('src/pdtime.r'); t0=Sys.time()
    if(args$dbvenn == 'venn') {
		venn_bed(b_path,out,fig_path,uni_list)
	} else if(args$dbvenn == 'summ') {
		summ_ann(b_path,anns,out)
	} else {
		paste0('[Error] There is no such function "',args$dbfilt,'" in db_filter: ',
            paste0(args$ldlink,collapse=', '),'\n') %>% cat
	}
	paste0(pdtime(t0,1),'\n') %>% cat
}