import subprocess
from string import Template

from AKRUP.funcbase import *

class RTemplate(object):
    """
    Creates a R script and runs it
    """

    def __init__(self, template, parameters):

        self.template = Template(template)
        self.parameters = parameters

    def run(self, clean=True):
        """
        Create a temporary file and run it
        """
        template = self.template
        parameters = self.parameters
        # write to a temporary R script
        fw = get_temname(mode='w')
        path = fw.name
        fw.write(template.safe_substitute(parameters))
        fw.close()
        # sh("Rscript %s" % path)
        try:
            p = subprocess.call(f'Rscript {path}', shell=True)
            if p != 0:
                print(f'Error!!! {p}')
        except Exception as e:
            print(e)

        if clean:
            os.remove(path)

class KsDistrubute:
    def __init__(self, options):

        self.args = dict(options)
        
    def run(self):
        self.rksdot(self.args)

    def rksdot(self, args):
        """
            run ks distrubute
        """

        dotks_template = """
            library(ggplot2)
            library(ggridges)

            mins=function(X, mu, sigma, a1){
                dnorm(X,mean=mu,sd=sigma)*a1
            }

            data<-read.csv("${ks_file}", header = T, encoding="UTF-8")  # T表示数据中的第一行是列名
            spec <- c()
            x <- c()
            y <- c()
            type <- c()
            color = c()
            ordername = c()
            len <- length(row.names(data))
            xl = seq(${area}, by = .001)
            num = 0
            for (i in 1:len){
                row_len = length(data[i,])
                row = na.omit(unlist(data[i,]))
                row_len = length(row)
                x <- c(x, xl)
                yl = mins(xl, data[i,3], data[i,4], data[i,2])
                type <- c(type, rep(paste('ty',num, sep=""), length(yl)))
                color <- c(color, data[i,5])
                ordername <- c(ordername, data[i,1])
                y <- c(y, yl)
                spec <- c(spec, rep(data[i,1], length(yl)))
                num <- num+1
            }
            data = data.frame(species = spec, x = x, y=y, type = type)

            lab = unique(data$type)
            data$species <-factor(data$species, levels = ordername)
            fig1 <- ggplot(data, aes(x = x, y = species, fill=type)) +
            geom_density_ridges(mapping = aes(height=y),inherit.aes=TRUE,stat = 'identity', scale = as.integer(${scale}), rel_min_height = 0.01, alpha=as.integer(${alpha}), color = 'NA')+
            theme_ridges(grid = FALSE)+
            scale_x_continuous(expand = c(0.01,0))+ # 扩展下横轴和纵轴
            scale_y_discrete(expand = c(0.01,0)) +
            theme(legend.position = "none", axis.text.x=element_text(hjust = 0.5, vjust = 0.5),axis.text.y=element_text(face="italic"),
                  axis.line.x.bottom=element_line(colour = "black", size=0.3), axis.line.y.left=element_line(colour = "black", size=0.3))+
            labs(x="ks",y="Kernel density of colinear block")+
            scale_fill_manual(values = color)
            ggsave("${save_fig}", width = as.integer(${width}), height = as.integer(${height}), dpi = as.integer(${dpi}))
        """

        rtemplate = RTemplate(dotks_template, args)
        rtemplate.run()
