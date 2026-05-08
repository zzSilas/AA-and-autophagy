library(igraph)
library(STRINGdb)
# 初始化 STRINGdb 对象（人类物种）
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")
options(timeout = 300)
mapped_data <- string_db$map(dif, "gene", removeUnmappedRows = TRUE)
hits <- mapped_data$STRING_id[1:200]  # 选择前200个基因
string_db$plot_network(hits)
colored_data <- string_db$add_diff_exp_color(subset(mapped_data, P.Value < 0.05), logFcColStr="logFC")
payload_id <- string_db$post_payload(colored_data$STRING_id, colors=colored_data$color)
string_db$plot_network(hits, payload_id=payload_id)

interactions <- string_db$get_interactions(hits)
setwd("E://斑秃")
write.csv(interactions, "ppi_data.csv", row.names = FALSE)
