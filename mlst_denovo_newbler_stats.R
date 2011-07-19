# parse newbler output files and report stats

# input are <454AllContigs.lengths> <454LargeContigs.lengths> <454Scaffolds.lengths> <outfile>
cat("-- reading arguments\n", sep = "")
Args <- commandArgs();

if (Args[5] == "NA") { 
   contig = read.table(file="454AllContigs.lengths", header=FALSE, sep="\t", as.is=TRUE)[,1]
   Lcontig = read.table(file="454LargeContigs.lengths", header=FALSE, sep="\t", as.is=TRUE)[,1]

   library(ggplot2)
   contig = read.table(file=Args[3], header=FALSE, sep="\t", as.is=TRUE)[,1]
   Lcontig = read.table(file=Args[4], header=FALSE, sep="\t", as.is=TRUE)[,1]

   # calc
   calcN50.f = function(con.v) {
      con.v = rev(sort(con.v))
      half = sum(con.v) / 2
      csum = cumsum(con.v)
      index = sum(csum <= half)
      con.v[index]
   }


   N50 = c(calcN50.f(contig), calcN50.f(Lcontig))
   AvgLen = round(c(mean(contig), mean(Lcontig)))
   MaxLen = c(max(contig), max(Lcontig))
   BasesCov = c(sum(contig), sum(Lcontig))

   # write out table with data

   df = data.frame(N50=N50, AvgLen = AvgLen, MaxLen = MaxLen, BasesCov=BasesCov)
   df.out = t(df)
   colnames(df.out) = c("AllContigs", "LargeContigs")

   write.table(df.out, Args[6], col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

   # Lcontigs
   Lcontig.df = data.frame(No=as.factor(1:length(Lcontig)), rev(sort(Lcontig)))
   colnames(Lcontig.df)[2] = "Lcontig"
   p = ggplot(Lcontig.df, aes(No, Lcontig)) + geom_histogram(fill="lightblue") + opts(axis.text.x=NULL, title="Large Contig (>500nt) lengths") + labs(x='', y='Length (nt)')
   p = p + geom_hline(aes(yintercept=N50[2]))
   pdf(file="contigLarge.lengths.pdf", width=7, height=7)
   print(p)
   dev.off()


} else {

   contig = read.table(file="454AllContigs.lengths", header=FALSE, sep="\t", as.is=TRUE)[,1]
   Lcontig = read.table(file="454LargeContigs.lengths", header=FALSE, sep="\t", as.is=TRUE)[,1]
   scaffolds = read.table(file="454Scaffolds.lengths", header=FALSE, sep="\t", as.is=TRUE)[,1]


   library(ggplot2)
   contig = read.table(file=Args[3], header=FALSE, sep="\t", as.is=TRUE)[,1]
   Lcontig = read.table(file=Args[4], header=FALSE, sep="\t", as.is=TRUE)[,1]
   scaffolds = read.table(file=Args[5], header=FALSE, sep="\t", as.is=TRUE)[,1]

   # calc
   calcN50.f = function(con.v) {
      con.v = rev(sort(con.v))
      half = sum(con.v) / 2
      csum = cumsum(con.v)
      index = sum(csum <= half)
      con.v[index]
   }


   N50 = c(calcN50.f(contig), calcN50.f(Lcontig), calcN50.f(scaffolds))
   AvgLen = round(c(mean(contig), mean(Lcontig), mean(scaffolds)))
   MaxLen = c(max(contig), max(Lcontig), max(scaffolds))
   BasesCov = c(sum(contig), sum(Lcontig), sum(scaffolds))

   # write out table with data

   df = data.frame(N50=N50, AvgLen = AvgLen, MaxLen = MaxLen, BasesCov=BasesCov)
   df.out = t(df)
   colnames(df.out) = c("AllContigs", "LargeContigs", "Scaffolds")

   write.table(df.out, Args[6], col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t")

   # plotting
   #contig.df = data.frame(No=as.factor(1:length(contig)), rev(sort(contig)))
   #colnames(contig.df)[2] = "contig"
   #p = ggplot(contig.df, aes(No, contig)) + geom_histogram(fill="lightblue") + opts(axis.text.x=NULL, title="All Contig (>100nt) lengths") + labs(x='', y='Length (nt)')
   #p = p + geom_hline(aes(yintercept=N50[1]))
   #pdf(file="contig.lengths.pdf", width=7, height=7)
   #print(p)
   #dev.off()

   # Lcontigs
   Lcontig.df = data.frame(No=as.factor(1:length(Lcontig)), rev(sort(Lcontig)))
   colnames(Lcontig.df)[2] = "Lcontig"
   p = ggplot(Lcontig.df, aes(No, Lcontig)) + geom_histogram(fill="lightblue") + opts(axis.text.x=NULL, title="Large Contig (>500nt) lengths") + labs(x='', y='Length (nt)')
   p = p + geom_hline(aes(yintercept=N50[2]))
   pdf(file="contigLarge.lengths.pdf", width=7, height=7)
   print(p)
   dev.off()

   # Scaffolds
   scaffolds.df = data.frame(No=as.factor(1:length(scaffolds)), rev(sort(scaffolds)))
   colnames(scaffolds.df)[2] = "scaffolds"
   p = ggplot(scaffolds.df, aes(No, scaffolds)) + geom_histogram(fill="lightblue") + opts(axis.text.x=NULL, title="Scaffolds (>500nt) lengths") + labs(x='', y='Length (nt)')
   p = p + geom_hline(aes(yintercept=N50[3]))
   pdf(file="scaffolds.lengths.pdf", width=7, height=7)
   print(p)
   dev.off()
}



