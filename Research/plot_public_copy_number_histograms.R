# Figure 3A: original histogram design, reading current saved AbsCNV results.
# No inference, re-fitting, re-segmentation, or CNV calling is performed here.
suppressWarnings(invisible(Sys.setlocale("LC_CTYPE", ".UTF-8")))
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
argv <- commandArgs(TRUE)
if(length(argv)!=3) stop("Usage: Rscript plot_public_copy_number_histograms.R META_DATA_CSV CN_PEAKS_CSV OUTPUT_DIR")
rd_path <- argv[1]; peak_path <- argv[2]; out <- argv[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(grepl("arial",systemfonts::match_fonts("Arial")$path[1],ignore.case=TRUE))
rd <- fread(rd_path)
peaks <- fread(peak_path)
samples <- c("MCF7","TAMR","FASR")
stopifnot(all(samples %in% names(rd)),all(c("sample","cn","mean") %in% names(peaks)))
stopifnot(setequal(unique(peaks$sample),samples))
stopifnot(all(rd[[1]] %in% paste0("chr",1:22)),all(rd[[3]]-rd[[2]]==1000))
histdata <- list(); histstats <- list()
for(s in samples) {
  # Original figure and Julia model use 50 consecutive retained rows, unchanged.
  smooth50 <- frollmean(rd[[s]],50,align="left")
  v <- smooth50[is.finite(smooth50)]
  stopifnot(all(v>=0 & v<=250))
  h <- hist(v,breaks=seq(0,250,2),plot=FALSE,include.lowest=TRUE)
  histdata[[s]] <- data.table(sample=s,x=h$mids,count=h$counts)
  histstats[[s]] <- data.table(sample=s,retained_bins=nrow(rd),valid_50bin_windows=length(v),excluded_windows=nrow(rd)-length(v))
}
hd <- rbindlist(histdata); hd[,sample:=factor(sample,levels=samples)]
pk <- copy(peaks)[,sample:=factor(sample,levels=samples)]
pk[,height:=mapply(function(ss,mu) max(hd[as.character(sample)==ss & abs(x-mu)<=6,count],0),as.character(sample),mean)]
pk <- merge(pk,hd[,.(maxcount=max(count)),by=sample],by="sample")
pk[,`:=`(label=sprintf("CN%d = %.1f",cn,mean),label_y=pmax(height*1.06,maxcount*.055))]
pk[as.character(sample)=="FASR" & cn==7,label_y:=maxcount*.18]
pk[as.character(sample)=="FASR" & cn==1,label_y:=maxcount*.33]
pk[,label_hjust:=fifelse(mean<40,0,fifelse(mean>220,1,.5))]
theme_set(theme_bw(base_size=8,base_family="Arial"))
theme_update(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=.2,color="#E8E8E8"),
 strip.background=element_rect(fill="#F0F2F4"),strip.text=element_text(size=8,face="bold"),
 plot.title=element_text(size=10,face="bold"),plot.subtitle=element_text(size=7.5),
 plot.caption=element_text(size=6.6,hjust=0),legend.position="bottom",legend.title=element_text(size=7),
 legend.text=element_text(size=7),plot.margin=margin(7,8,6,7))
ph <- ggplot(hd,aes(x,count))+geom_col(width=2,fill="#9ECAE1",color="#6BAED6",linewidth=.12)+
 geom_segment(data=pk,aes(x=mean,xend=mean,y=0,yend=height),inherit.aes=FALSE,color="#B71920",linewidth=.45)+
 geom_text(data=pk,aes(x=mean,y=label_y,label=label,hjust=label_hjust),inherit.aes=FALSE,color="#A50F15",size=2.4,family="Arial",fontface="bold",vjust=0)+
 facet_wrap(~sample,nrow=1,scales="free_y")+scale_x_continuous(limits=c(0,270),breaks=seq(0,250,50),expand=c(0,0))+
 scale_y_continuous(expand=expansion(mult=c(0,.2)),labels=scales::label_number(scale=1e-3,suffix="k"))+
 labs(x="Normalized read depth",y="Frequency",title=NULL,subtitle=NULL,caption=NULL)
stem <- file.path(out,"Figure3A_public_CN_histograms")
ggsave(paste0(stem,".pdf"),ph,width=170,height=62,units="mm",device=cairo_pdf,family="Arial",bg="white")
ggsave(paste0(stem,".svg"),ph,width=170,height=62,units="mm",device=svglite::svglite,bg="white")
ragg::agg_png(paste0(stem,".png"),width=170,height=62,units="mm",res=600,background="white");print(ph);dev.off()
fwrite(hd,file.path(out,"histogram_bins.csv"),bom=TRUE)
fwrite(rbindlist(histstats),file.path(out,"histogram_window_audit.csv"),bom=TRUE)
fwrite(pk,file.path(out,"plotted_peak_labels.csv"),bom=TRUE)
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
cat("Figure 3A updated from current saved RD and GMM parameters.\n")
