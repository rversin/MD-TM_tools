# Get data
dt_tmp = read.table("contact.dat", header = TRUE)

contact = data.frame(atom1 = dt_tmp$atom1, atom2 = dt_tmp$atom2, count = dt_tmp$count)

heatmap = ggplot(data = contact, aes(x=atom1, y=atom2, fill=count)) +
  geom_tile(color = "white") +
  ylab("TM2") +
  xlab("TM1") +
  scale_fill_gradientn(colors = c("#CCCCCC", "#FFFF33", "#FFCC33", "#FF3300", "#FF3399", "#660033"),
                       limit = c(0.00001, 1),   
                       breaks = c(0.00001, 0.2, 0.4, 0.6, 0.8, 1),
                       space = "Lab", name = "") +
  theme(legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, hjust = 1, face = "bold"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 14, hjust = 1, face = "bold"),
        strip.text.x = element_text(size = 20, face = "bold")) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 20,
                               title.position = "top", title.hjust = 0.5)) +
  coord_fixed()

ggsave("contact-map.png", heatmap, width = 5, height = 5, limitsize = FALSE)