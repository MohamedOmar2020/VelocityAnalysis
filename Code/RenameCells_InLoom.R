library(loomR)

ldata1 <- connect("~/Desktop/Looms/MAD1.loom", mode = "r+", skip.validate = T)
ldata2 <- connect("~/Desktop/Looms/MAD2.loom", mode = "r+", skip.validate = T)
ldata3 <- connect("~/Desktop/Looms/MAD3.loom", mode = "r+", skip.validate = T)
ldata4 <- connect("~/Desktop/Looms/MAD4.loom", mode = "r+", skip.validate = T)
ldata5 <- connect("~/Desktop/Looms/MAD5.loom", mode = "r+", skip.validate = T)
ldata6 <- connect("~/Desktop/Looms/MAD6.loom", mode = "r+", skip.validate = T)

ldata1[["col_attrs/obs_names"]][] <- paste(ldata1[["col_attrs/obs_names"]][], "-1", sep = "")
ldata2[["col_attrs/obs_names"]][] <- paste(ldata2[["col_attrs/obs_names"]][], "-2", sep = "")
ldata3[["col_attrs/obs_names"]][] <- paste(ldata3[["col_attrs/obs_names"]][], "-3", sep = "")
ldata4[["col_attrs/obs_names"]][] <- paste(ldata4[["col_attrs/obs_names"]][], "-4", sep = "")
ldata5[["col_attrs/obs_names"]][] <- paste(ldata5[["col_attrs/obs_names"]][], "-5", sep = "")
ldata6[["col_attrs/obs_names"]][] <- paste(ldata6[["col_attrs/obs_names"]][], "-6", sep = "")






## Close the loom connection
ldata1$close_all()
ldata2$close_all()
ldata3$close_all()
ldata4$close_all()
ldata5$close_all()
ldata6$close_all()


###
