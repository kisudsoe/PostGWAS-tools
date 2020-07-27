# 190807 Time result print to cat.
# 191107 Revert the result as return 
pdtime = function(time,type=1) {
    t=Sys.time()
    d=difftime(t,time,unit="sec")
    if(d<60) d=paste(round(d,1),"sec")
    else if(d>=60 && d<3600) d=paste(round(d/60,1),"min")
    else if(d>=3600) d=paste(round(d/3600,1),"hr")
    if(type==1) out = paste("Job done:",t,"for",d)
    else if(type==2) out = paste("Job process:",d)
    return(out)
}