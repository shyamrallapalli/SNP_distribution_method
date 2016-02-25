#encoding: utf-8

require 'rinruby'

class Plot
  def self.qqplot(experimental, dir, title, ylabel, xlabel, nametag)
    myr = RinRuby.new(:echo=>false)
    myr.assign 'experimental', experimental
    myr.assign 'title', title
    myr.assign 'xlabel', xlabel
    myr.assign 'ylabel', ylabel
    myr.assign 'dir', dir
    myr.assign 'nametag', nametag
    myr.eval 'qqline2 <- function(x, y, probs = c(0.25, 0.75), qtype = 7, ...){
      stopifnot(length(probs) == 2)
      x2 <- quantile(x, probs, names=FALSE, type=qtype, na.rm = TRUE)
      y2 <- quantile(y, probs, names=FALSE, type=qtype, na.rm = TRUE)
      slope <- diff(y2)/diff(x2)
      int <- y2[1L] - slope*x2[1L]
      abline(int, slope, ...)
    }

    leg_r2 <- function(k){
      legend(x = "topleft", bty = "n",
             legend = substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                                 list(a = format(coef(k)[1], digits = 2),
                                      b = format(coef(k)[2], digits = 2),
                                      r2 = format(summary(k)$r.squared, digits = 3))))
    }

    pdf(paste(dir, "/", nametag, "_qqplot_exp_hyp.pdf", sep=""), width=4, height=4)
    par(cex.axis=0.5, cex.lab=0.8, cex.main=1, mar=c(2.5,2,1,0.2), oma=c(0,0,0,0), mgp=c(1, 0.3, 0))
    options(scipen = 10)
    x <- experimental
    y <- rnorm(length(x), mean(x), sd(x))
    # df <- data.frame(x, y)
    V = qqplot(x, y, main=title, ylab=ylabel, xlab =xlabel)
    l <- qqline2(x, y, col = 6)
    # fg <- data.frame(V$x, V$y)
    k <- lm(V$y ~ V$x)
    len <- leg_r2(k)
    dev.off()'
    myr.quit
  end

  def self.densities(hm, ht, ratio, dir)
    myr = RinRuby.new(:echo=>false)
    myr.assign 'hm', hm
    myr.assign 'ht', ht
    myr.assign 'ratio', ratio
    myr.assign 'dir', dir
    myr.eval 'pdf(paste(dir, "/experimental_densities.pdf", sep=""), width=7, height=3)
    par(cex.axis=0.5, cex.lab=0.8, cex.main=1, mar=c(2.5,2,1,0.2), oma=c(0,0,0,0), mgp=c(1, 0.3, 0))
    options(scipen = 10)
    d1 <- density(hm, adjust = 1, kernel = c("gaussian"))
    d2 <- density(ht, adjust = 1, kernel = c("gaussian"))
    d3 <- density(ratio, adjust =1 , kernel = c("gaussian"))
    # limit the plot on x-axis based on variant position spread
    length <- range(hm, ht, ratio)
    plot(range(d1$x, d2$x, d3$x), range(d1$y, d2$y, d3$y), type = "n",
      xlim =c(length[1], length[2]), xlab = "variant position", ylab = "densities")
    lines(d1, col = "magenta2") ## Homozygous
    lines(d2, col = "royalblue2", lty=2) ## Heterozygous
    lines(d3, col = "gray46", lwd =3) ## Homozygous/Heterozygous ratio
    axis(3, at=c(length[1], length[1]+diff(length)/2, length[2]), labels=c(1, diff(length)/2, diff(length)) )
    legend("topright",col=c("magenta2", "royalblue2", "grey46"),lwd=1,lty=1:2,
      legend=c("Homozygous SNP density","Heterozygous SNP density", "Hom/Het ratio"), bty="n")
    dev.off()'
    myr.quit
  end

  def self.get_ylim(array, genome_length)
    myr = RinRuby.new(:echo=>false)
    myr.assign 'array', array
    myr.assign 'genome_length', genome_length
    myr.eval 'plot((1:512)*(genome_length/512), density(array)$y)'
    ylim = myr.pull 'par("yaxp")[2] + par("yaxp")[2]/par("yaxp")[3]'
    myr.quit
    ylim
  end

  def self.comparison(real_ratios, experimental, length, dir, ylim, original_pos, outcome_pos)
    myr = RinRuby.new(:echo=>false)
    myr.assign 'experimental', experimental
    myr.assign 'real_ratios', real_ratios
    myr.assign 'length', length
    myr.assign 'dir', dir
    myr.assign 'ylim', ylim
    myr.assign 'original_pos', original_pos
    myr.assign 'outcome_pos', outcome_pos

    myr.eval 'd1 <-density(experimental, adjust = 1, kernel = c("gaussian"))
    d2 <-density(real_ratios, adjust = 1, kernel = c("gaussian"))
    pdf(paste(dir, "/compare_outcome.pdf", sep=""), width=7, height=4)
    par(cex.axis=0.6, cex.lab=1, cex.main=1.2, mar=c(2.5,1.2,1,0.2), oma=c(0,0,0,0), mgp=c(1, 0.3, 0))
    layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(1,3))
    #options(scipen = 10)
    plot(d1$x, d1$y, col = "slategray4", lwd =2, type = "l", main="sdm ratios", xlab="", ylab="")
    plot(d1$x, d1$y, type = "n", main = "variant arrangment", xlim =c(0,length), xlab="", ylab="", yaxt="n")
    abline(v=outcome_pos, col = "blue", lty=2)
    plot(d2$x, d2$y, col = "steelblue3", lwd =2, type = "l", main="original ratios", xlab="", ylab="")
    plot(d2$x, d2$y, type = "n", main = "original var positions", xlim =c(0,length), xlab="genome positions", ylab="", yaxt="n")
    abline(v=original_pos, col = "red", lty=2)
    dev.off()'
    myr.quit
  end

end
