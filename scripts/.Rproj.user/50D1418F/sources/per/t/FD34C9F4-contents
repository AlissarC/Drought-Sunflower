
install.packages("gridExtra")
library(gridExtra)
library(grid)


############################## Fig 1 a and b ##################################################

FC_Fig <- read.csv('../input/FC.csv')
names(FC_Fig)

FC_Fig$plant = as.character(FC_Fig$plant)

B2= subset(FC_Fig, ID_trt=="B2")
B5= subset(FC_Fig, ID_trt=="B5")
B6= subset(FC_Fig, ID_trt=="B6")
B3= subset(FC_Fig, ID_trt=="B3")

SPD= subset(FC_Fig, drought=="SPD")
FDD= subset(FC_Fig, drought=="FDD")


Fig1a <- ggplot(data = B5, aes(x = DAD, y = FC, color = plant)) +
  geom_point(size = 3) +   # Add points with shapes for nitrogen levels
  geom_line(aes(group = plant), size = 1)+
  theme(legend.position = "right",
        plot.title = element_text(size = rel(5), hjust = 1, vjust = 1),  
        plot.title.position = "panel",  
        legend.title = element_text(size = rel(2.2)),
        legend.text = element_text(size = rel(2.2)),
        plot.tag = element_text(size = rel(4)),
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_text(size=rel(4), colour = 'black'),
        axis.text.y=element_text(size=rel(4), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  ylab(expression(italic('FC')[''] *' (g water g soil' ^ '-1' * ')')) +
  xlab(expression (italic('Day After Drought')))+
  labs(title = "HN-SPD")  # Add title here
Fig1a



Fig1b <- ggplot(data = B2, aes(x = DAD, y = FC, color = plant)) +
  geom_point(size = 3) +   # Add points with shapes for nitrogen levels
  geom_line(aes(group = plant), size = 1)+
  theme(legend.position = "right",
        plot.title = element_text(size = rel(5), hjust = 1, vjust = 1),  
        plot.title.position = "panel",  
        legend.title = element_text(size = rel(2.2)),
        legend.text = element_text(size = rel(2.2)),
        plot.tag = element_text(size = rel(4)),
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_text(size=rel(4), colour = 'black'),
        axis.text.y=element_text(size=rel(4), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  ylab(expression(italic('FC')[''] *' (g water g soil' ^ '-1' * ')')) +
  xlab(expression (italic('Day After Drought')))+
  labs(title = "LN-SPD")  # Add title here
Fig1b


Fig1c <- ggplot(data = B6, aes(x = DAD, y = FC, color = plant)) +
  geom_point(size = 3) +   # Add points with shapes for nitrogen levels
  geom_line(aes(group = plant), size = 1)+
  theme(legend.position = "right",
        plot.title = element_text(size = rel(5), hjust = 1, vjust = 1),  
        plot.title.position = "panel",  
        legend.title = element_text(size = rel(2.2)),
        legend.text = element_text(size = rel(2.2)),
        plot.tag = element_text(size = rel(4)),
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_text(size=rel(4), colour = 'black'),
        axis.text.y=element_text(size=rel(4), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  ylab(expression(italic('FC')[''] *' (g water g soil' ^ '-1' * ')')) +
  xlab(expression (italic('Day After Drought')))+
  labs(title = "HN-FDD")  # Add title here
Fig1c


Fig1d <- ggplot(data = B3, aes(x = DAD, y = FC, color = plant)) +
  geom_point(size = 3) +   # Add points with shapes for nitrogen levels
  geom_line(aes(group = plant), size = 1)+
  theme(legend.position = "right",
        plot.title = element_text(size = rel(5), hjust = 1, vjust = 1),  
        plot.title.position = "panel",  
        legend.title = element_text(size = rel(2.2)),
        legend.text = element_text(size = rel(2.2)),
        plot.tag = element_text(size = rel(4)),
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_text(size=rel(4), colour = 'black'),
        axis.text.y=element_text(size=rel(4), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  ylab(expression(italic('FC')[''] *' (g water g soil' ^ '-1' * ')')) +
  xlab(expression (italic('Day After Drought')))+
  labs(title = "LN-FDD")  # Add title here
Fig1d

Fig1a.g <- ggplotGrob(Fig1a)
Fig1b.g <- ggplotGrob(Fig1b)
Fig1c.g <- ggplotGrob(Fig1c)
Fig1d.g <- ggplotGrob(Fig1d)
SPD <- cbind (Fig1a.g,Fig1b.g, size="max" )
FDD <- cbind (Fig1c.g,Fig1d.g, size="max" )
Fig1 <- rbind(SPD,FDD,size = 'max')

jpeg(filename = "Fig1.png", 
     width = 16, height = 14, units = 'in', res = 1400)
grid.newpage()
grid.draw(Fig1)

grid.text("(a)", x = 0.1, y = 0.98, just = "left", gp = gpar(fontsize = 24, fontface = "bold"))
grid.text("(b)", x = 0.60, y = 0.98, just = "left", gp = gpar(fontsize = 24, fontface = "bold"))
grid.text("(c)", x = 0.1, y = 0.48, just = "left", gp = gpar(fontsize = 24, fontface = "bold"))
grid.text("(d)", x = 0.60, y = 0.48, just = "left", gp = gpar(fontsize = 24, fontface = "bold"))

dev.off()
############################ panel1 Fig4 ##########################################
gs_plot_t4.g <- ggplotGrob(gs_plot_t4)
Chi_plot_t4.g <- ggplotGrob(Chi_plot_t4)
betaleaf_plot_t4.g <- ggplotGrob(betaleaf_plot_t4)
panel1_Fig4 <- rbind (gs_plot_t4.g,Chi_plot_t4.g, betaleaf_plot_t4.g,size="max")

Nmass_plot_t4.g <- ggplotGrob(Nmass_plot_t4)
LMA_plot_t4.g <- ggplotGrob(LMA_plot_t4)
Narea_plot_t4.g <- ggplotGrob(Narea_plot_t4)
panel2_Fig4 <- rbind (Nmass_plot_t4.g,LMA_plot_t4.g, Narea_plot_t4.g,size="max")


Fig4 <- cbind(panel1_Fig4,panel2_Fig4,size = 'max')
                                                                   
jpeg(filename = "Fig4.png", 
     width = 14, height = 16, units = 'in', res = 1200)
grid.newpage()
grid.draw(Fig4)

grid.text("(a)", x = 0.03, y = 0.99, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(b)", x = 0.03, y = 0.67, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(c)", x = 0.025, y = 0.33, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

grid.text("(d)", x = 0.93, y = 0.99, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(e)", x = 0.93, y = 0.67, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(f)", x = 0.93, y = 0.33, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

# Close the JPEG device
dev.off()

############################ Fig 5  ##########################################
Vcmax_plot_t4.g <- ggplotGrob(Vcmax_plot_t4)
Anet_plot_t4.g <- ggplotGrob(Anet_plot_t4)
panel1_Fig5 <- rbind (Vcmax_plot_t4.g,Anet_plot_t4.g, size="max" )

PNUE_plot_t4.g <- ggplotGrob(PNUE_plot_t4)
WUE_plot_t4.g <- ggplotGrob(WUE_plot_t4)
panel2_Fig5 <- rbind (WUE_plot_t4.g,PNUE_plot_t4.g, size="max" )

Fig5 <- cbind(panel1_Fig5,panel2_Fig5,size = 'max')

jpeg(filename = "Fig5.png", 
     width = 16, height = 14, units = 'in', res = 1200)
grid.newpage()
grid.draw(Fig5)

grid.text("(a)", x = 0.013, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(b)", x = 0.013, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

grid.text("(c)", x = 0.93, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(d)", x = 0.93, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

# Close the JPEG device
dev.off()

####################### Fig6 ##############################################################
PSA_plot_t4.g <- ggplotGrob(PSA_plot_t4)
biomass_plot_t4.g <- ggplotGrob(biomass_plot_t4)
panel1_Fig6 <- cbind (PSA_plot_t4.g,biomass_plot_t4.g, size="max")

N_cost_plot_t4_plot.g <- ggplotGrob(N_cost_plot_t4)
transp_cost_plot_t4.g <- ggplotGrob(transp_cost_plot_t4)
panel2_Fig6 <- cbind (N_cost_plot_t4_plot.g,transp_cost_plot_t4.g, size="max")

Fig6 <- rbind(panel1_Fig6,panel2_Fig6,size = 'max')

jpeg(filename = "Fig6.png", 
     width = 16, height = 14, units = 'in', res = 1400)
grid.newpage()
grid.draw(Fig6)

grid.text("(a)", x = 0.02, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(b)", x = 0.93, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(c)", x = 0.02, y = 0.48, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(d)", x = 0.93, y = 0.48, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

# Close the JPEG device
dev.off()


#######################################################################################
############################ Fig S3  ##########################################
Centered_gs_plot.g <- ggplotGrob(Centered_gs_plot)
Centered_Chi_plot.g <- ggplotGrob(Centered_Chi_plot)
Centered_Betaleaf_plot.g <- ggplotGrob(Centered_Betaleaf_plot)

panel1_FigS3 <- rbind (Centered_gs_plot.g, Centered_Chi_plot.g, Centered_Betaleaf_plot.g,size="max")

Centred_Nmass_plot.g <- ggplotGrob(Centred_Nmass_plot)
Centered_LMA_plot.g <- ggplotGrob(Centered_LMA_plot)
Centred_Narea_plot.g <- ggplotGrob(Centered_Narea_plot)
panel2_FigS3 <- rbind (Centred_Nmass_plot.g, Centered_LMA_plot.g, Centred_Narea_plot.g,size="max")


FigS3 <- cbind(panel1_FigS3,panel2_FigS3,size = 'max')

jpeg(filename = "FigS3.png", 
     width = 14, height = 16, units = 'in', res = 1200)
grid.newpage()
grid.draw(FigS3)

grid.text("(a)", x = 0.02, y = 0.99, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(b)", x = 0.03, y = 0.67, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(c)", x = 0.03, y = 0.33, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

grid.text("(d)", x = 0.93, y = 0.99, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(e)", x = 0.93, y = 0.67, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(f)", x = 0.93, y = 0.33, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

# Close the JPEG device
dev.off()

############################ Fig S4 ##########################################
Centred_Vcmax_plot.g <- ggplotGrob(Centred_Vcmax_plot)
Centered_Anet_plot.g <- ggplotGrob(Centered_Anet_plot)
panel1_FigS4 <- rbind (Centred_Vcmax_plot.g,Centered_Anet_plot.g, size="max")

Centered_PNUE_plot.g <- ggplotGrob(Centered_PNUE_plot)
Centered_LWUE_plot.g <- ggplotGrob(Centered_LWUE_plot)
panel2_FigS4 <- rbind (Centered_LWUE_plot.g,Centered_PNUE_plot.g, size="max")

FigS4 <- cbind(panel1_FigS4,panel2_FigS4)

jpeg(filename = "FigS4.png", 
     width = 16, height = 14, units = 'in', res = 1200)
grid.newpage()
grid.draw(FigS4)

grid.text("(a)", x = 0.013, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(b)", x = 0.013, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

grid.text("(c)", x = 0.93, y = 0.98, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("(d)", x = 0.93, y = 0.50, just = "left", gp = gpar(fontsize = 14, fontface = "bold"))

# Close the JPEG device
dev.off()


##################### Figure S1 #########################################
fig_S1_subset = subset(stat_sunflower, drought %in% c("SPD", "FDD"))
FigS1 <- ggplot(data= fig_S1_subset, aes(x=drought, y= Sdt_DD, color = N_level))+
  geom_boxplot() +
  geom_jitter(aes(shape = drought), size = 3) +
  scale_shape_manual(values = c("FDD" = 17, "SPD" = 15)) +
  theme(legend.position = "right",
        plot.title = element_text(size = rel(2.2)),
        legend.title = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        plot.tag = element_text(size = rel(2)),
        axis.title.y=element_text(size=rel(2.2), colour = 'black'),
        axis.title.x=element_text(size=rel(2.2), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  
  ylab(expression (italic('Drought Severity Index')))+
  xlab(expression (italic('Drought scenarios')))
FigS1

ggsave("FigS1.jpeg", plot = FigS1, 
       width = 30, height = 20, units = "cm")

############################ Temp and PAR ##########################################
Temp_plot.g <- ggplotGrob(Temp)
PAR_plot.g <- ggplotGrob(PAR)

Greenhouse <- rbind (Temp_plot.g,PAR_plot.g, size="max" )

jpeg(filename = "Greenhouse.jpeg", 
     width = 14, height = 16, units = 'in', res = 600)
grid.newpage()
grid.draw(Greenhouse)
grid.text("(a)", x = 0.03, y = 0.99, just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
grid.text("(b)", x = 0.03, y = 0.50, just = "left", gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()



