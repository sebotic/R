#script retrieves melanoma incidence rates from Statistics Austria and generates a publication quality figure with ggplot

require(ggplot2)
require(XML)

url <- "http://www.statistik.at/web_de/statistiken/gesundheit/krebserkrankungen/haut/021736.html"


tables <- readHTMLTable(url)

mel.incidence <- as.data.frame(tables[[1]], stringsAsFactors=F)

colnames(mel.incidence) <- c("year", "total", "male", "female",
                             "total (age standardized)", "male (age standardized)", "female (age standardized)",
                            "total (cumulative)", "male (cumulative)", "female (cumulative)")

#quick and dirty conversion from german to english locale
for(i in 1:ncol(mel.incidence)){
  mel.incidence[, i] <- gsub("\\.", "", mel.incidence[, i])
  mel.incidence[, i] <- as.numeric(gsub(",", ".", mel.incidence[, i]))
}

p <- ggplot(data = mel.incidence, aes(x = year)) + 
  geom_line(aes(y = male, linetype="male")) +
  geom_line(aes(y = female, linetype="female")) +
  geom_line(aes(y = total, linetype="total")) + 
  xlab("Year") +
  ylab("Patients per year") +
  ggtitle("Melanoma Incidence from 1983 - 2011 in Austria (Data: Statistics Austria)") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  coord_cartesian(xlim=c(1983, 2011)) + 
  scale_x_continuous(breaks=seq(1983, 2011, 2)) +
  scale_colour_manual(values=c("black", "black", "black")) +  
  scale_linetype_manual(values=c("dashed", "dotted", "solid"))
  
p 
