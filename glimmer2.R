
#change the number for the glimmer predict
plotGlimmer = function(file = '17.glimmer.predict') {
t = read.table(file, header = F, skip = 1)
c = data.frame(size=abs(t[,2]-t[,3]))
ggplot(c, aes(size))+geom_histogram(binwidth = 1000)+ggtitle(file)
}
#change the png name
ggsave('17.png')
plotGlimmer()

