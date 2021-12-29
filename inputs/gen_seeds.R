seed_seq <- sample(0:1248, 1000)
write.table(seed_seq, file = "seed_samples.csv", row.names = FALSE, col.names = FALSE)
